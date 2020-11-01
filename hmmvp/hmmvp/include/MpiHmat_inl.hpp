/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

// Routines to parallelize Hmat::Mvp().

#include <stdio.h>
#include <algorithm>
#include "util/include/CodeAnalysis.hpp"
#include "util/include/Mpi.hpp"
#include "hmmvp/include/Hmat.hpp"
#include "hmmvp/include/HmatIo.hpp"

namespace hmmvp {

template<typename T> Blint MpiHmat<T>::GetM () const { return hm->GetM(); }
template<typename T> Blint MpiHmat<T>::GetN () const { return hm->GetN(); }
template<typename T> Blint MpiHmat<T>::GetNnz () const { return nnz; }
template<typename T> const Blint* MpiHmat<T>::GetQ () const
{ return hm->GetQ(); }
template<typename T> const Blint* MpiHmat<T>::GetP () const
{ return hm->GetP(); }

namespace {
struct HmatData { UInt m, n, realp, nb, nnz; };
struct BlockData { Blint i, nnz, r0, c0; };

template<bool sort_by_col>
bool CmpBlockDataCR (const BlockData& bd1, const BlockData& bd2) {
  if (sort_by_col)
    return (bd1.c0 == bd2.c0) ? bd1.r0 < bd2.r0 : bd1.c0 < bd2.c0;
  else
    return (bd1.r0 == bd2.r0) ? bd1.c0 < bd2.c0 : bd1.r0 < bd2.r0;
}

bool CmpBlockDataNnz (const BlockData& bd1, const BlockData& bd2) {
  return bd1.nnz > bd2.nnz;
}
      
void GetBlocks (const string& filename, HmatData* h, vector<BlockData>* bds) {
  double tol;
  h->nnz = 0;
  FILE* fid = fopen(filename.c_str(), "rb");

  { Blint bm, bn, brealp, bnb;
    hmmvp::ReadHmatHeader(fid, bm, bn, brealp, bnb, tol);
    if (bm <= 0 || bn <= 0 || !(brealp == 1 || brealp == 2) || bnb <= 0)
      throw FileException("Problem in reading header.");
    h->m = (UInt) bm; h->n = (UInt) bn; h->realp = (UInt) brealp;
    h->nb = (UInt) bnb; }

  bds->resize(h->nb);
  for (UInt i = 0; i < h->nb; i++) {
    Blint bm, bn, brank;
    BlockData& bd = (*bds)[i];
    bd.i = i;
    hmmvp::ReadHmatBlockInfo(fid, h->realp, bd.r0, bd.c0, bm, bn, brank);
    if (brank == std::min(bm, bn)) bd.nnz = bm*bn;
    else bd.nnz = (bm + bn)*brank;
    h->nnz += bd.nnz;
  }
  fclose(fid);        
}

Blint SetupMpi (UInt nproc, const string& filename,
                vector<vector<Blint> >& block_idxs)
  throw (FileException, Exception)
{
  HmatData hd;
  vector<BlockData> bds_ns;
  GetBlocks(filename, &hd, &bds_ns);
  assert(hd.nb == bds_ns.size());
  if (hd.nb < nproc)
    throw Exception("Reduce number of MPI processes. Currently too many "
                    "for this problem.");

  // Sort the blocks by decreasing nnz.
  std::sort(bds_ns.begin(), bds_ns.end(), CmpBlockDataNnz);

  // Allocate blocks among tasks draft style: from 1 to nproc, then nproc to
  // 1, and repeat.
  vector<UInt> seps(nproc+1);
  vector<BlockData> bds(bds_ns.size());
  for (size_t ip = 0, ib = 0; ip < nproc; ip++) {
    seps[ip] = ib;
    const size_t ind1 = 2*(nproc - ip) - 1, ind2 = 2*ip + 1;
    for (size_t in = ip, ibi = 0; in < bds.size() && ib < bds.size();
         ib++, ibi++) {
      bds[ib] = bds_ns[in];
      // This carries out the draft ordering.
      in += (ibi % 2 == 0) ? ind1 : ind2;
    }
  }
  seps[nproc] = bds.size();

  block_idxs.clear();
  block_idxs.resize(nproc);
  for (UInt i = 0; i < nproc; i++) {
    // Sort each task's blocks by increasing row. This should encourage
    // cache coherence in y space in each task.
    std::sort(&bds[seps[i]], &bds[seps[i+1]], CmpBlockDataCR<false>);
    // Create block_idxs array.
    block_idxs[i].clear();
    block_idxs[i].resize(seps[i+1] - seps[i]);
    for (UInt j = seps[i], k = 0; j < seps[i+1]; j++, k++)
      block_idxs[i][k] = bds[j].i;
  }

  return hd.nnz;
}

static const int tag0 = 100;
}

template<typename T>
MpiHmat<T>::MpiHmat (const string& filename, UInt ncol, UInt nthreads)
  throw (FileException, Exception)
  : workr(0), am_root(false), do_perms(true)
{
  int pid = mpi::Pid();

  am_root = mpi::AmRoot();
  vector<Blint> my_block_idxs;
  if (am_root) {
    int nproc = mpi::GetNproc();
    // Determine how to distribute blocks among processes.
    vector<vector<Blint> > block_idxs;
    try {
      nnz = SetupMpi(nproc, filename, block_idxs);
      mpi_IsTrue(true);
    } catch (const Exception& e) {
      mpi_IsTrue(false);
      throw e;
    }
    // Assign everyone else their blocks.
    for (int i = 1; i < nproc; i++) {
      Blint nb = block_idxs[i].size();
      mpi::Send(&nb, 1, i, tag0);
      mpi::Send(&block_idxs[i][0], nb, i, tag0);
    }
    // Assign myself my blocks.
    my_block_idxs = block_idxs[0];
  } else {
    // Wait until I receive my block indices.
    Blint nb;
    if (!mpi_IsTrue()) throw Exception();
    mpi::Recv(&nb, 1, 0, tag0);
    my_block_idxs.resize(nb);
    mpi::Recv(&my_block_idxs[0], nb, 0, tag0);
  }
  
  // Now read in just my blocks.
  try {
    hm = hmmvp::NewHmat(filename, ncol, nthreads, &my_block_idxs);
  } catch (const hmmvp::FileException& fe) {
    throw FileException(fe.GetMsg());
  }
  hm->TurnOffPermute();

  UInt nw = ncol*std::max(hm->GetM(), hm->GetN());
  workr.resize(nw);
}

template<typename T> MpiHmat<T>::~MpiHmat () { hmmvp::DeleteHmat(hm); }

template<typename T>
void MpiHmat<T>::Mvp (const T* x, T* y, UInt ncol) const
  throw (Exception)
{
  if (do_perms) {
    // For all but root, x can be NULL. (Not advertised in interface.)
    if (am_root) hm->ApplyQ(x, &workr[0], ncol);
    // For now, send Q*x to all procs.
    mpi::Bcast(&workr[0], ncol*hm->GetN());
    // Do my part.
    hm->Mvp(&workr[0], y, ncol);
    // Sum all the y's.
    mpi::Reduce(y, &workr[0], ncol*hm->GetM(), MPI_SUM);
    if (am_root) {
      // Permute back to user's basis.
      hm->ApplyPt(&workr[0], y, ncol);
    }
  } else {
    // Do my part. (NB: x must be valid in every task, unlike the other case.)
    hm->Mvp(x, &workr[0], ncol);
    // Sum all the y's.
    mpi::Reduce(&workr[0], y, (int) ncol*hm->GetM(), MPI_SUM);
  }
}

}
