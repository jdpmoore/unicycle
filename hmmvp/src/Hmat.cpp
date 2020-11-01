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

// TODO:
// - handle ncol at Block level rather than Hmat level
//   x done for no-(rs,cs) case using gemm rather than gemv
// - better-scaling parallelization
//   x seem to have something good for Pre/PostMult and Mvp for no-(rs,cs) case
// - Transpose MVP
//   x done for no-(rs,cs) case
// - change memsets for e2b case. provide optional &rs, &cs args to Pre/PostMult.
// - make sure arithmetic is >= ~50% of the runtime for various cases.

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include "util/include/Util.hpp"
#include "util/include/OpenMP.hpp"
#include "util/include/LinAlg.hpp"
#include "hmmvp/include/HmatIo.hpp"
#include "hmmvp/include/Hmat.hpp"
using namespace util;

// For debugging purposes: This makes arithmetic proceed in the same order for
// various test cases.
#define SORT_BLOCKS 0

//#define TV
#ifdef TV
#include <sys/time.h>

static double difftime (struct timeval& t1, struct timeval& t2) {
  static const double us = 1.0e6;
  return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
}
#endif

namespace hmmvp {
#define DELETE(a) if (a) delete[] (a);

typedef vector< vector<Blint> > VecVecBlint;
template<typename T> struct Block;
template<typename Real> class TypedHmat;

// For element-to-block maps.
struct E2BMap {
  VecVecBlint is[2];
  UInt nb;
  // flag is used to indicate whether a block has been included.
  typedef UInt FlagType;
  FlagType flag;
  FlagType* use_block;

  E2BMap() : flag(0), use_block(0) {}
  ~E2BMap() { DELETE(use_block); }

  template<typename Real>
  void Init(const TypedHmat<Real>& hm);
  void InitFlag();
  UInt GetUsedBlocks(UInt cdim, const Blint* cs, UInt ncs, Blint* blks);
  UInt GetUsedBlocks(UInt cdim, const Blint* cs, UInt ncs,
                     UInt rdim, const Blint* rs, UInt nrs,
                     Blint* blks);
};

// Implement single and double precision operations.
template<typename Real> class TypedHmat : public Hmat {
  friend class HmatConstructor;

public:
  TypedHmat(const string& filename, UInt ncol = 1, UInt max_nthreads = 1,
            // Optionally read in only blocks block_idxs[0]:nb_idxs[1]-1.
            const vector<Blint>* block_idxs = 0)
    throw (FileException);
  ~TypedHmat();

  virtual UInt GetScalarSize () const { return sizeof(Real); }
  virtual UInt SetThreads(UInt nthreads);
  virtual void ReorganizeMemory();

  virtual double NormFrobenius2() const;
  virtual double NormOne() const;

  virtual UInt GetNnz(const vector<Blint>& rs, const vector<Blint>& cs);

  virtual void Mvp (const float *x, float *y, UInt ncol)
    throw (Exception) { Mvp<float>(x, y, ncol); }
  virtual void Mvp (const double *x, double *y, UInt ncol)
    throw (Exception) { Mvp<double>(x, y, ncol); }
  virtual void MvpT (const float *x, float *y, UInt ncol)
    throw (Exception) { MvpT<float>(x, y, ncol); }
  virtual void MvpT (const double *x, double *y, UInt ncol)
    throw (Exception) { MvpT<double>(x, y, ncol); }

  virtual void Mvp (const float* x, float* y, UInt ncol,
                    const vector<Blint>* rs,
                    const vector<Blint>* cs)
    throw (Exception)
  { Mvp<float>(x, y, ncol, rs, cs); }
  virtual void Mvp (const double* x, double* y, UInt ncol,
                    const vector<Blint>* rs,
                    const vector<Blint>* cs)
    throw (Exception)
  { Mvp<double>(x, y, ncol, rs, cs); }

  virtual void ApplyQ (const float* x, float* xq, UInt ncol) const {
    ApplyQ<float>(x, xq, ncol); }
  virtual void ApplyQ (const double* x, double* xq, UInt ncol) const {
    ApplyQ<double>(x, xq, ncol); }

  virtual void ApplyPt (const float* yp, float* y, UInt ncol) const {
    ApplyPt<float>(yp, y, ncol); }
  virtual void ApplyPt (const double* yp, double* y, UInt ncol) const {
    ApplyPt<double>(yp, y, ncol); }

  virtual const Blint* GetQ () const { return q; }
  virtual const Blint* GetP () const { return p; }

  virtual void Extract (const vector<Blint>& rs, const vector<Blint>& cs,
                        float* es) {
    Extract<float>(rs, cs, es); }
  virtual void Extract(const vector<Blint>& rs, const vector<Blint>& cs,
                       double* es) {
    Extract<double>(rs, cs, es); }

  virtual void FullBlocksIJS (vector<Blint>& I, vector<Blint>& J,
                              vector<float>& S, float cutoff) const {
    FullBlocksIJS<float>(I, J, S, cutoff); }
  virtual void FullBlocksIJS (vector<Blint>& I, vector<Blint>& J,
                              vector<double>& S, double cutoff) const {
    FullBlocksIJS<double>(I, J, S, cutoff); }

  Blint GetNbrBlocks () const { return nb; }
  const vector<Block<Real> >& GetBlocks () const { return blocks; }

private:
  UInt K;
  UInt nb; // nbr blocks
  UInt max_nthreads, nthreads;
  UInt ncol;

  vector<Block<Real> > blocks;
  Real* blocks_mem;
  Real* work;
  Real* xq;
  Real* ys;
  Blint* p; // row permutation (base-0 indexing)
  Blint* q; // column permutation

  // For B(rs,cs) MVPs
  Blint *ip, *iq; // inverse permutations
  Blint *rs, *cs, *rs_idxs, *cs_idxs;
  Blint* blks;
  UInt kblks;
  // For Extract(). Inverse perms of subset of rows and cols.
  Blint *irs, *ics;
  E2BMap e2bm;

  // Indices in row and column space for each thread
  vector<Blint> thr_row_seps, thr_col_seps;
  vector<Blint> thr_block_seps;

  bool transp_mode;

private:
  // For HmatConstructor.
  TypedHmat(UInt nb, const vector<Blint>& p, const vector<Blint>& q,
            UInt ncol = 1, UInt max_nthreads = 1);

  template<typename T>
  void Mvp(const T *x, T *y, UInt ncol) throw (Exception);
  template<typename T>
  void RMvp(const T *x, T *y, UInt ncol, const vector<Blint>& rs)
    throw (Exception);
  template<typename T>
  void CMvp(const T *x, T *y, UInt ncol, const vector<Blint>& cs)
    throw (Exception);
  template<typename T>
  void Mvp(const T *x, T *y, UInt ncol, const vector<Blint>& rs,
           const vector<Blint>& cs) throw (Exception);
  template<typename T>
  void Mvp(const T *x, T *y, UInt ncol, const vector<Blint>* rs,
           const vector<Blint>* cs) throw (Exception);
  template<typename T>
  void MvpT(const T *x, T *y, UInt ncol) throw (Exception);

  template<typename T> void ApplyQ(const T* x, T* xq, UInt ncol) const;
  template<typename T> void ApplyPt(const T* yp, T* y, UInt ncol) const;

  template<typename T>
  void Extract(const vector<Blint>& rs, const vector<Blint>& cs, T* es);

  template<typename T>
  void FullBlocksIJS(vector<Blint>& I, vector<Blint>& J, vector<T>& S,
                     T cutoff = (T) -1.0) const;

  void ReadHmatFile(FILE* fid, const vector<Blint>* block_idxs = 0)
    throw (FileException);
  void SetupArrays();

  template<typename T> void PreMult(const T* x, UInt nthreads, UInt ncol);
  template<typename T> void PostMult(UInt nthreads, UInt ncol, T* y);
  template<typename T>
  void PreMult(const T* x, UInt nthreads, UInt ncol, const vector<Blint>& cs);
  template<typename T>
  void PreMult(const T* x, UInt nthreads, UInt ncol, const vector<Blint>& rs,
               const vector<Blint>& cs);
  template<typename T> void PostMult(UInt nthreads, UInt ncol, T* y,
                                     const vector<Blint>& rs);
  void InvPerms();
  bool RIdxLimits(
    const Block<Real>& b, UInt nrs, Blint& rlo, Blint& rhi, Blint& nbrs);
  void RIdxLimitsNoTest(
    const Block<Real>& b, UInt nrs, Blint& rlo, Blint& rhi, Blint& nbrs);
  bool CIdxLimits(
    const Block<Real>& b, UInt ncs, Blint& clo, Blint& chi, Blint& nbcs);
  void CIdxLimitsNoTest(
    const Block<Real>& b, UInt ncs, Blint& clo, Blint& chi, Blint& nbcs);
  bool IdxLimits(
    const Block<Real>& b, UInt nrs, UInt ncs, Blint& rlo, Blint& rhi,
    Blint& clo, Blint& chi, Blint& nbrs, Blint& nbcs);
  void IdxLimitsNoTest(
    const Block<Real>& b, UInt nrs, UInt ncs, Blint& rlo, Blint& rhi,
    Blint& clo, Blint& chi, Blint& nbrs, Blint& nbcs);
  void SetIdxMapsE2B(UInt nrs, UInt ncs);
  void RGetUsedBlocks(UInt nrs);
  void CGetUsedBlocks(UInt ncs);
  void GetUsedBlocks(UInt nrs, UInt ncs);
  UInt GetUsedBlocksE2B(UInt nrs, UInt ncs);
  void RMvpBlock(UInt ncol, const Block<Real>& b, UInt nrs);
  void CMvpBlock(UInt ncol, const Block<Real>& b, UInt ncs);
  void MvpBlock(UInt ncol, const Block<Real>& b, UInt nrs, UInt ncs);
  template<typename T>
  void ExtractBlock(const Block<Real>& b, Blint nrs, Blint ncs, T* es);

  void SetupParallel();
  void SwapForTranspose();
};

template<typename Real> template<typename T> inline void TypedHmat<Real>::
Mvp (const T *x, T *y, UInt ncol, const vector<Blint>* rs,
     const vector<Blint>* cs)
  throw (Exception)
{
  if (rs == NULL) {
    if (cs == NULL) Mvp(x, y, ncol);
    else CMvp(x, y, ncol, *cs);
  } else {
    if (cs == NULL) RMvp(x, y, ncol, *rs);
    else Mvp(x, y, ncol, *rs, *cs);
  }
}

};

namespace hmmvp {
template<typename Real>
struct Block {
  // r0, c0 are 0-based
  Blint m, n, r0, c0;
  // rank if U,Vt are used
  Blint rank;
  // [U (with s mult'ed in) and V^T] or [dense block B]
  Real *U, *Vt, *B;
  bool am_responsible_for_mem;

  Block();

  int Read(FILE* fid) throw (FileException);
  int Set(Blint r0, Blint c0, const Matrix<Real>& B);
  int Set(Blint r0, Blint c0, const Matrix<Real>& U, const Matrix<Real>& V);
  ~Block();
  void Reallocate();
  size_t UseThisMemory(Real* mem);
  void Delete();

  void Mvp(const bool transp, const Real* x, Real* y, Real* work) const;
  void Mvp(const bool transp, const Real* x, Real* y, Real* work,
           UInt ldy, UInt ldx, UInt nrhs) const;
  void RMvp(const Real* x, Real* y, Real* work, const Blint* rs, Blint nrs)
    const;
  void CMvp(const Real* x, Real* y, Real* work, const Blint* cs, Blint ncs)
    const;
  void Mvp(const Real* x, Real* y, Real* work, const Blint* rs, Blint nrs,
           const Blint* cs, Blint ncs) const;

  double NormFrobenius2() const;
  void AccumColNorm1s(Real* work, vector<double>& col_norm1s) const;
  double MaxMag() const;

  UInt Nnz() const;
  UInt Nnz(Blint nrs, Blint ncs) const;
};
  
template<typename Real>
Block<Real>::Block () : U(0), Vt(0), B(0), am_responsible_for_mem(true) {}

// Read the blocks file to construct a block. Return nnz.
template<typename Real>
int Block<Real>::Read (FILE *fid) throw (FileException) {
  int nnz = -1;
  try {
    ReadHmatBlock(fid, r0, c0, m, n, rank, B, U, Vt);
    nnz = B ? m*n : (m + n)*rank;
  } catch(FileException& e) {
    Delete();
    throw e;
  }
  return nnz;
}

template<typename Real>
int Block<Real>::Set (Blint r0, Blint c0, const Matrix<Real>& B) {
  this->r0 = r0;
  this->c0 = c0;
  m = B.Size(1);
  n = B.Size(2);
  rank = std::min(m, n);
  this->B = new Real[m*n];
  memcpy(this->B, B.GetPtr(), m*n*sizeof(Real));
  U = 0;
  Vt = 0;
  return m*n;
}

template<typename Real> int Block<Real>::
Set (Blint r0, Blint c0, const Matrix<Real>& U, const Matrix<Real>& V) {
  assert(U.Size(2) == V.Size(2));
  this->r0 = r0;
  this->c0 = c0;
  m = U.Size(1);
  n = V.Size(1);
  rank = U.Size(2);
  B = 0;
  this->U = new Real[m*rank];
  memcpy(this->U, U.GetPtr(), m*rank*sizeof(Real));
  this->Vt = new Real[rank*n];
  Matrix<Real> mVt(rank, n, this->Vt);
  Transpose(V, mVt);
  return rank*(m + n);
}

template<typename Real>
inline void ReallocateArray (Real*& v, const size_t n, bool& arfm) {
  if (!v) return;
  Real* vnew = new Real[n];
  memcpy(vnew, v, n*sizeof(Real));
  if (arfm) delete[] v;
  v = vnew;
  arfm = true;
}

template<typename Real> void Block<Real>::Reallocate () {
  ReallocateArray(B, m*n, am_responsible_for_mem);
  ReallocateArray(U, m*rank, am_responsible_for_mem);
  ReallocateArray(Vt, n*rank, am_responsible_for_mem);
}

template<typename Real>
size_t UseThisMemoryForArray (Real* mem, Real*& v, const size_t n, bool& arfm) {
  if (!v) return 0;
  memcpy(mem, v, n*sizeof(Real));
  if (arfm) delete[] v;
  v = mem;
  arfm = false;
  return n;
}

template<typename Real> size_t Block<Real>::UseThisMemory (Real* mem) {
  size_t k;
  if ((k = UseThisMemoryForArray(mem, B, m*n, am_responsible_for_mem)))
    return k;

  size_t k1, k2;
  bool arfm = am_responsible_for_mem;
  k1 = UseThisMemoryForArray(mem, Vt, n*rank, arfm);
  mem += k1;
  k2 = UseThisMemoryForArray(mem, U, m*rank, am_responsible_for_mem);
  return k1 + k2;
}

template<typename Real> Block<Real>::~Block () { Delete(); }

template<typename Real>
void Block<Real>::Delete () {
  if (!am_responsible_for_mem) return;
  DELETE(U); DELETE(Vt); DELETE(B);
}

template<typename Real> UInt Block<Real>::Nnz () const {
  if (B) return m*n;
  else return (m + n)*rank;
}

template<typename Real> UInt Block<Real>::Nnz (Blint nrs, Blint ncs) const {
  if (B) return nrs*ncs;
  else return (nrs + ncs)*rank;
}

// gemm-based MVP. Not currently used, so it may have fallen into disrepair,
// especially the transp block.
template<typename Real> inline void Block<Real>::
Mvp (const bool transp, const Real* x, Real* y, Real* work,
     UInt ldy, UInt ldx, UInt nrhs) const {
  const Block& b = *this;
  if (transp) {
    if (b.B) {
      // y(c0+(0:n-1),:) += B'*x(r0+(0:m-1),:)
      gemm('t', 'n', b.n, nrhs, b.m,
           (Real) 1.0, b.B, b.m, x + b.r0, ldx,
           (Real) 1.0, y + b.c0, ldy);
    } else {
      if (b.rank == 0) return;
      // work(1:rank,:) += U'*x(r0+(0:m-1),:)
      gemm('t', 'n', b.rank, nrhs, b.m,
           (Real) 1.0, b.U, b.m, x + b.r0, ldx,
           (Real) 1.0, work, b.rank);
      // y(c0+(0:n-1),:) = Vt'*work(1:rank,:)
      gemm('t', 'n', b.n, nrhs, b.rank,
           (Real) 1.0, b.Vt, b.rank, work, b.rank,
           (Real) 0.0, y + b.c0, ldy);
    }
  } else {
    if (b.B) {
      // y(r0+(0:m-1),:) += B*x(c0+(0:n-1),:)
      gemm('n', 'n', b.m, nrhs, b.n,
           (Real) 1.0, b.B, b.m, x + b.c0, ldx,
           (Real) 1.0, y + b.r0, ldy);
    } else {
      if (b.rank == 0) return;
      // work(1:rank,:) = Vt*x(c0+(0:n-1),:)
      gemm('n', 'n', b.rank, nrhs, b.n,
           (Real) 1.0, b.Vt, b.rank, x + b.c0, ldx,
           (Real) 0.0, work, b.rank);
      // y(r0+(0:m-1),:) += U*work(1:rank,:)
      gemm('n', 'n', b.m, nrhs, b.rank,
           (Real) 1.0, b.U, b.m, work, b.rank,
           (Real) 1.0, y + b.r0, ldy);
    }
  }
}

template<typename Real> inline void Block<Real>::
Mvp (const bool transp, const Real* x, Real* y, Real* work) const {
  const Block& b = *this;
  if (transp) {
    if (b.B) {
      // y(r0+(0:m-1)) += B*x(c0+(0:n-1))
      gemv('t', b.m, b.n,
           (Real) 1.0, b.B, b.m, x + b.r0, 1,
           (Real) 1.0, y + b.c0, 1);
    } else {
      if (b.rank == 0) return;
      // work(1:rank) += U'*x(r0+(0:m-1))
      gemv('t', b.m, b.rank,
           (Real) 1.0, b.U, b.m, x + b.r0, 1,
           (Real) 0.0, work, 1);
      // y(c0+(0:n-1)) = Vt'*work(1:rank)
      gemv('t', b.rank, b.n,
           (Real) 1.0, b.Vt, b.rank, work, 1,
           (Real) 1.0, y + b.c0, 1);
    }
  } else {
    if (b.B) {
      // y(r0+(0:m-1)) += B*x(c0+(0:n-1))
      gemv('n', b.m, b.n,
           (Real) 1.0, b.B, b.m, x + b.c0, 1,
           (Real) 1.0, y + b.r0, 1);
    } else {
      if (b.rank == 0) return;
      // work(1:rank) = Vt*x(c0+(0:n-1))
      gemv('n', b.rank, b.n,
           (Real) 1.0, b.Vt, b.rank, x + b.c0, 1,
           (Real) 0.0, work, 1);
      // y(r0+(0:m-1)) += U*work(1:rank)
      gemv('n', b.m, b.rank,
           (Real) 1.0, b.U, b.m, work, 1,
           (Real) 1.0, y + b.r0, 1);
    }
  }
}

template<typename Real> inline void Block<Real>::
RMvp (const Real* x, Real* y, Real* work, const Blint* rs, Blint nrs) const {
  const Block& b = *this;
  Real* yp = y + b.r0;
  const Real* xp = x + b.c0;
  if (b.B) {
    // I've decided not to use 'dot'. The arithmetic comes out differently for
    // some reason. In principle it's not a problem since it's at machine
    // precision, but in practice I like the arithmetic to be the same for
    // debugging purposes.
#if 0
    for (Blint i = 0; i < nrs; i++) {
      Blint r = rs[i] - b.r0;
      yp[r] += dot(b.n, b.B + r, b.m, xp, 1);
    }
#else
    const Real* B = b.B;
    for (Blint c = 0; c < b.n; c++) {
      Real xc = xp[c];
      for (Blint i = 0; i < nrs; i++) {
        Blint r = rs[i] - b.r0;
        yp[r] += B[r]*xc;
      }
      B += b.m;
    }
#endif
  } else {
    if (b.rank == 0) return;
    // work(1:rank) = Vt*x(c0+(0:n-1))
    gemv('n', b.rank, b.n,
         (Real) 1.0, b.Vt, b.rank, xp, 1,
         (Real) 0.0, work, 1);
    // y(r0+rs) += U(rs,:)*work(1:rank)
#if 0
    for (Blint i = 0; i < nrs; i++) {
      Blint r = rs[i] - b.r0;
      yp[r] += dot(b.rank, b.U + r, b.m, work, 1);
    }
#else
    const Real* U = b.U;
    for (Blint j = 0; j < b.rank; j++) {
      Real wj = work[j];
      for (Blint i = 0; i < nrs; i++) {
        Blint r = rs[i] - b.r0;
        yp[r] += U[r]*wj;
      }
      U += b.m;
    }
#endif
  }
}

template<typename Real> inline void Block<Real>::
CMvp (const Real* x, Real* y, Real* work, const Blint* cs, Blint ncs) const {
  const Block& b = *this;
  Real* yp = y + b.r0;
  const Real* xp = x + b.c0;
  if (b.B) {
    for (Blint i = 0; i < ncs; i++) {
      Blint c = cs[i] - b.c0;
      axpy(b.m, xp[c], b.B + b.m*c, 1, yp, 1);
    }
  } else {
    if (b.rank == 0) return;
    // work(1:rank) = Vt(:,cs)*x(c0+cs)
    memset(work, 0, b.rank*sizeof(Real));
    for (Blint i = 0; i < ncs; i++) {
      Blint c = cs[i] - b.c0;
      axpy(b.rank, xp[c], b.Vt + b.rank*c, 1, work, 1);
    }
    // y(r0+(0:m-1)) += U*work(1:rank)
    gemv('n', b.m, b.rank,
         (Real) 1.0, b.U, b.m, work, 1,
         (Real) 1.0, yp, 1);
  }
}

template<typename Real> inline void Block<Real>::
Mvp (const Real* x, Real* y, Real* work, const Blint* rs, Blint nrs,
     const Blint* cs, Blint ncs) const {
  const Block& b = *this;
  const Real* xp = x + b.c0;
  Real* yp = y + b.r0;
  if (b.B) {
    // y(r0+rs) = B(rs,cs)*x(cs)
    for (Blint j = 0; j < ncs; j++) {
      Blint c = cs[j] - b.c0;
      const Real* B = b.B + c*b.m;
      Real xc = xp[c];
      for (Blint i = 0; i < nrs; i++) {
        Blint r = rs[i] - b.r0;
        yp[r] += B[r]*xc;
      }
    }
  } else {
    if (b.rank == 0) return;
    // work(1:rank) = Vt(:,cs)*x(c0+cs)
    UInt m = b.rank;
    memset(work, 0, m*sizeof(Real));
    for (Blint j = 0; j < ncs; j++) {
      Blint c = cs[j] - b.c0;
      const Real* Vt = b.Vt + c*m;
      Real xc = xp[c];
      for (UInt i = 0; i < m; i++)
        work[i] += Vt[i]*xc;
    }    
    // y(r0+rs) += U(rs,:)*work(1:rank)
    const Real* U = b.U;
    for (UInt j = 0; j < m; j++) {
      Real wj = work[j];
      for (Blint i = 0; i < nrs; i++) {
        Blint r = rs[i] - b.r0;
        yp[r] += U[r]*wj;
      }
      U += b.m;
    }
  }
}

template<typename Real> inline double Block<Real>::NormFrobenius2 () const {
  const Block& b = *this;
  double nf2 = 0.0;
  if (b.B) {
    Blint mn = b.m*b.n;
    for (Blint i = 0; i < mn; i++)
      nf2 += B[i]*B[i];
  } else {
    if (b.rank == 0) return nf2;
    // Could decide to do sum(sum(b.U'.*((b.V'*b.V)*b.U'))) based on
    // b.m,n. Not now.
    char transa, transb;
    blas_int m, n, k, lda, ldb, ldc;
    Real alpha, beta;
    // nf2 = sum(sum(b.V'.*((b.U'*b.U)*b.V')))
    Blint r2 = b.rank*b.rank, rn = b.rank*b.n;
    vector<Real> work(r2 + rn);
    // work = U' U
    transa = 't';
    lda = b.m;
    m = n = b.rank;
    k = b.m;
    transb = 'n';
    ldb = b.m;
    ldc = b.rank;
    alpha = 1.0;
    beta = 0.0;
    gemm(transa, transb, m, n, k, alpha, b.U, lda, b.U, ldb,
         beta, &work[0], ldc);
    // work = work V'
    transa = 'n';
    lda = b.rank;
    n = b.n;
    k = b.rank;
    ldb = b.rank;
    Real* wp = &work[0] + r2;
    gemm(transa, transb, m, n, k, alpha, &work[0], lda, b.Vt, ldb,
         beta, wp, ldc);
    // sum(V.*work)
    for (Blint i = 0; i < rn; i++)
      nf2 += wp[i]*b.Vt[i];
  }
  return nf2;
}

template<typename Real> inline void Block<Real>::
AccumColNorm1s (Real* work, vector<double>& col_norm1s) const {
  const Block& b = *this;
  if (b.B) {
    Blint os = 0;
    for (Blint j = 0; j < b.n; j++) {
      for (Blint i = 0; i < b.m; i++)
        col_norm1s[b.c0 + j] += fabs(b.B[os + i]);
      os += b.m;
    }
  } else {
    Blint os = 0;
    for (Blint j = 0; j < b.n; j++) {
      // work = U*Vt(:,j)
      gemv('n', b.m, b.rank, (Real) 1.0, b.U, b.m,
           b.Vt + os, 1,
           (Real) 0.0, work, 1);
      for (Blint i = 0; i < b.m; i++)
        col_norm1s[b.c0 + j] += (double) fabs(work[i]);
      os += b.rank;
    }
  }
}

template<typename Real> double MaxMag (const Real* a, Blint N) {
  double max_mag = 0.0;
  for (Blint i = 0; i < N; i++)
    if (fabs(a[i]) > max_mag) max_mag = fabs(a[i]);
  return max_mag;
}

template<typename Real> inline double Block<Real>::MaxMag () const {
  const Block& b = *this;
  if (b.B)
    return hmmvp::MaxMag(b.B, b.m*b.n);
  else
    return hmmvp::MaxMag(b.U, b.m*b.rank)*hmmvp::MaxMag(b.Vt, b.rank*b.n);
}

Hmat* NewHmat (const string& filename, UInt ncol, UInt imax_nthreads,
               const vector<Blint>* block_idxs)
  throw (FileException)
{
  Blint m, n, realp, nb;
  double tol;
  HmatInfo(filename, m, n, realp, nb, tol);
  if (realp == 1)
    return new TypedHmat<float>(filename, ncol, imax_nthreads, block_idxs);
  else if (realp == 2)
    return new TypedHmat<double>(filename, ncol, imax_nthreads, block_idxs);
  else
    throw FileException(string("NewHmat: realp is not 1 or 2 in ") + filename);
}

void DeleteHmat (Hmat* hm) { delete hm; }

template<typename T> static double ReadBlockAndGetNormFrobenius2 (FILE* fid) {
  Block<T> b;
  b.Read(fid);
  return b.NormFrobenius2();
}

double HmatNormFrobenius2 (const string& filename) throw (FileException) {
  FILE* fid = fopen(filename.c_str(), "rb");
  if (!fid) throw FileException(string("Hmat: can't read file ") + filename);
  Blint m, n, realp, nb;
  double tol;
  ReadHmatHeader(fid, m, n, realp, nb, tol);
  double nf2 = 0.0;
  for (Blint i = 0; i < nb; i++)
    nf2 += (realp == 1) ?
      ReadBlockAndGetNormFrobenius2<float> (fid) :
      ReadBlockAndGetNormFrobenius2<double>(fid);
  fclose(fid);
  return nf2;
}

Hmat::Hmat ()
  : m(0), n(0), nnz(0), do_perms(true), save_state(false), state_saved(false),
    do_e2b_map(false)
{}

Hmat::~Hmat () {}

UInt Hmat::GetM () const { return m; }
UInt Hmat::GetN () const { return n; }
UInt Hmat::GetNnz () const { return nnz; }
void Hmat::TurnOnPermute ()  { do_perms = true;  state_saved = false; }
void Hmat::TurnOffPermute () { do_perms = false; state_saved = false; }
void Hmat::SaveState ()    { save_state = true; }
void Hmat::ReleaseState () { save_state = false; state_saved = false; }
void Hmat::TurnOnElementToBlockMap ()  { do_e2b_map = true;  }
void Hmat::TurnOffElementToBlockMap () { do_e2b_map = false; }

HmatConstructor::HmatConstructor (
  UInt realp, UInt nb, const vector<Blint>& p, const vector<Blint>& q,
  UInt ncol, UInt max_nthreads) throw (Exception)
  : _realp(realp), _hmf(NULL), _hmd(NULL), _called_GetHmat(false)
{
  if (realp == 1)
    _hmf = new TypedHmat<float>(nb, p, q, ncol, max_nthreads);
  else if (realp == 2)
    _hmd =  new TypedHmat<double>(nb, p, q, ncol, max_nthreads);
  else
    throw FileException("realp must be 1 or 2");
}
HmatConstructor::~HmatConstructor () {
  if (!_called_GetHmat) {
    if (_hmf) delete _hmf;
    else if (_hmd) delete _hmd;
  }
}
void HmatConstructor::
SetBlockB (UInt bi, UInt r0, UInt c0, const Matrix<double>& B) {
  assert(_realp == 2);
  assert(bi < _hmd->nb);
  _hmd->nnz += _hmd->blocks[bi].Set(r0, c0, B);
}
void HmatConstructor::
SetBlockB (UInt bi, UInt r0, UInt c0, const Matrix<float>& B) {
  assert(_realp == 1);
  assert(bi < _hmf->nb);
  _hmf->nnz += _hmf->blocks[bi].Set(r0, c0, B);
}
void HmatConstructor::
SetBlockUV (UInt bi, UInt r0, UInt c0, const Matrix<double>& U,
            const Matrix<double>& V) {
  assert(_realp == 2);
  assert(bi < _hmd->nb);
  _hmd->nnz += _hmd->blocks[bi].Set(r0, c0, U, V);
}
void HmatConstructor::
SetBlockUV (UInt bi, UInt r0, UInt c0, const Matrix<float>& U,
            const Matrix<float>& V) {
  assert(_realp == 1);
  assert(bi < _hmf->nb);
  _hmf->nnz += _hmf->blocks[bi].Set(r0, c0, U, V);
}
Hmat* HmatConstructor::GetHmat () {
  if (!_called_GetHmat) {
    if (_realp == 1) {
      _hmf->SetThreads(_hmf->max_nthreads);
      _hmf->SetupArrays();
    } else {
      _hmd->SetThreads(_hmd->max_nthreads);
      _hmd->SetupArrays();
    }
  }
  _called_GetHmat = true;
  if (_realp == 1) return _hmf;
  else return _hmd;
}

template<typename Real> TypedHmat<Real>::TypedHmat (
  const string& filename, UInt ncol, UInt imax_nthreads,
  const vector<Blint>* block_idxs)
  throw (FileException)
  : ncol(ncol), blocks_mem(0), work(0), xq(0), ys(0), p(0), q(0), ip(0), iq(0),
    rs(0), cs(0), rs_idxs(0), cs_idxs(0), blks(0), irs(0), ics(0),
    transp_mode(false)
{
  if (this->ncol == 0) this->ncol = 1;

  FILE* fid = fopen(filename.c_str(), "rb");
  if (!fid) throw FileException(string("Hmat: can't read file ") + filename);
  ReadHmatFile(fid, block_idxs);
  fclose(fid);

  max_nthreads = imax_nthreads;
  SetThreads(imax_nthreads);

  SetupArrays();
}

template<typename Real> TypedHmat<Real>::
TypedHmat (UInt nb, const vector<Blint>& p, const vector<Blint>& q, UInt ncol,
           UInt max_nthreads)
  : nb(nb), ncol(ncol), blocks_mem(0), work(0), xq(0), ys(0), p(0), q(0), ip(0),
    iq(0), rs(0), cs(0), rs_idxs(0), cs_idxs(0), blks(0), irs(0), ics(0),
    transp_mode(false)
{
  if (this->ncol == 0) this->ncol = 1;
  this->max_nthreads = max_nthreads;

  m = p.size();
  n = q.size();
  this->p = new Blint[m];
  this->q = new Blint[n];
  for (UInt i = 0; i < m; i++) this->p[i] = p[i];
  for (UInt i = 0; i < n; i++) this->q[i] = q[i];

  blocks.resize(nb);
}

template<typename Real> TypedHmat<Real>::~TypedHmat () {
  DELETE(blocks_mem);
  DELETE(work); DELETE(xq); DELETE(ys);
  DELETE(p); DELETE(q);
  DELETE(ip); DELETE(iq);
  DELETE(rs); DELETE(cs); DELETE(rs_idxs); DELETE(cs_idxs);
  DELETE(blks);
  DELETE(irs); DELETE(ics);
}

template<typename Real> void TypedHmat<Real>::SetupArrays () {
  // Set up other stuff.
  K = std::min(m, n);
  work = new Real[nthreads*ncol*K];
  xq = new Real[nthreads*ncol*n];
  ys = new Real[nthreads*ncol*m];
}

namespace {
struct BlockIdx {
  Blint idx, i;
  bool operator< (const BlockIdx& bi) const { return idx < bi.idx; }
};
}

template<typename Real> void TypedHmat<Real>::
ReadHmatFile (FILE* fid, const vector<Blint>* block_idxs)
  throw (FileException)
{
  // Read the header.
  { vector<FileBlint> vp, vq;
    double tol;
    Blint realp, bm, bn, bnb;
    ReadHmatHeader(fid, bm, bn, realp, bnb, tol, vp, vq);
    if (bm <= 0 || bn <= 0 || bnb <= 0)
      throw FileException("Hmat: Size is <= 0.");
    m = (UInt) bm; n = (UInt) bn; nb = (UInt) bnb;
    if ((realp == 1 && sizeof(Real) != sizeof(float)) ||
        (realp == 2 && sizeof(Real) != sizeof(double)) ||
        realp < 1 || realp > 2)
      throw FileException("Hmat: Wrong real precision.");
    p = new Blint[m];
    q = new Blint[n];
    for (UInt i = 0; i < m; i++) p[i] = (Blint) vp[i];
    for (UInt i = 0; i < n; i++) q[i] = (Blint) vq[i]; }

  // Read the blocks.
  nnz = 0;
  if (block_idxs) {
    // Read in a list of blocks. This is for MpiHmat.
    nb = block_idxs->size();
    // Sort block indices in increasing order. We want the permutation.
    vector<BlockIdx> sort_block_idxs(nb);
    for (UInt i = 0; i < nb; i++) {
      sort_block_idxs[i].idx = (*block_idxs)[i];
      sort_block_idxs[i].i = i;
    }
    std::sort(sort_block_idxs.begin(), sort_block_idxs.end());
    // Read in the desired blocks.
    blocks.resize(nb);
    for (Blint i = 0, k = 0; i <= sort_block_idxs[nb-1].idx; i++)
      if (i == sort_block_idxs[k].idx) {
        // We want this block.
        Block<Real>& b = blocks[sort_block_idxs[k].i];
        nnz += b.Read(fid);
        k++;
      } else {
        // Skip.
        Blint br0, bc0, bm, bn, brank;
        ReadHmatBlockInfo<Real>(fid, br0, bc0, bm, bn, brank);
      }
  } else {
    // Read in all the blocks.
    blocks.resize(nb);
    for (UInt i = 0; i < nb; i++) {
      Block<Real>& b = blocks[i];
      nnz += b.Read(fid);
    }
  }
}

template<typename Real> UInt TypedHmat<Real>::SetThreads (UInt nthreads) {
  if (nthreads == 0) nthreads = 1;
  if (nthreads > max_nthreads) nthreads = max_nthreads;
  if (nthreads > 1) {
    omp_set_num_threads(nthreads);
    this->nthreads = omp_get_max_threads();
    assert(this->nthreads <= nthreads);
  } else this->nthreads = nthreads;
  SetupParallel();
  return this->nthreads;
}

inline void SetThreadIndexSeps (vector<Blint>& seps, UInt nthreads, UInt n) {
  seps.clear();
  seps.resize(nthreads + 1);
  UInt nc = n / nthreads + 1;
  seps[0] = 0;
  for (UInt ti = 0; ti < nthreads - 1; ti++)
    seps[ti+1] = (ti + 1)*nc;
  seps[nthreads] = n;
}

#if 0
// I think this method is not scaling well to many (> 4, certainly > 8) OpenMP
// threads.
template<typename Real> void TypedHmat<Real>::SetupParallel () {
  // Assign blocks to threads. I tried sorting them in some way, but it's far
  // more important to have blocks contiguous in memory. The strategy here is
  // to split blocks into groups having roughly equal total nnz.
  thr_block_seps.clear();
  thr_block_seps.resize(nthreads + 1);
  thr_block_seps[0] = 0;
  assert(nthreads >= 1);
  for (UInt i = 0, cnnz = 0, tid = 0, nnzpt = nnz / nthreads + 1;
       i < nb && tid < nthreads - 1; i++) {
    cnnz += blocks[i].Nnz();
    if (cnnz > nnzpt) {
      thr_block_seps[tid+1] = i;
      cnnz = 0;
      tid++;
    }
  }
  thr_block_seps[nthreads] = nb;

  // Divide permuted rows and cols into contiguous group. Cols belonging to
  // thread tid are
  //     thr_cols_seps[tid]:thr_col_seps[tid+1]-1
  // and similarly for rows.
  SetThreadIndexSeps(thr_row_seps, nthreads, m);
  SetThreadIndexSeps(thr_col_seps, nthreads, n);
}
#else
// New in version 1.1.
struct BlockData { Blint i, nnz, r0; };
static bool CmpBlockDataNnz (const BlockData& bd1, const BlockData& bd2) {
  return bd1.nnz > bd2.nnz;
}

template<typename Real> void TypedHmat<Real>::SetupParallel () {
  assert(nthreads >= 1);

  vector<BlockData> bds(blocks.size());
  { // Sort the blocks by decreasing nnz.
    vector<BlockData> bds_ns(blocks.size());
    for (size_t i = 0; i < bds_ns.size(); i++) {
      const Block<Real>& b = blocks[i];
      bds_ns[i].i = i;
      bds_ns[i].nnz = b.B ? b.m*b.n : b.rank*(b.m + b.n);
    }
    std::sort(bds_ns.begin(), bds_ns.end(), CmpBlockDataNnz);

    // Allocate blocks among tasks draft style: from 1 to nthreads, then nthreads
    // to 1, and repeat.
    thr_block_seps.clear();
    thr_block_seps.resize(nthreads + 1);
    for (size_t ip = 0, ib = 0; ip < nthreads; ip++) {
      thr_block_seps[ip] = ib;
      const size_t ind1 = 2*(nthreads - ip) - 1, ind2 = 2*ip + 1;
      for (size_t in = ip, ibi = 0; in < bds.size() && ib < bds.size();
           ib++, ibi++) {
        bds[ib] = bds_ns[in];
        // This carries out the draft ordering.
        in += (ibi % 2 == 0) ? ind1 : ind2;
      }
    }
    thr_block_seps[nthreads] = bds.size();
  }

#ifndef NDEBUG
  { vector<char> is_used(bds.size(), 0);
    for (size_t i = 0; i < bds.size(); i++) {
      assert(!is_used[i]);
      is_used[i] = true;
    }
    for (size_t i = 0; i < bds.size(); i++) assert(is_used[i]); }
#endif

  { // Now rebuild the blocks array. Block's copy constructor is the default and
    // so shallow. Hence I'm just moving pointers around.
    vector<Block<Real> > old_blocks(blocks);
    for (size_t i = 0; i < blocks.size(); i++)
      blocks[i] = old_blocks[bds[i].i];
    // But watch out for the destructor.
    for (size_t i = 0; i < old_blocks.size(); i++) {
      Block<Real>& b = old_blocks[i];
      b.B = b.U = b.Vt = 0;
    }
  }

  // Divide permuted rows and cols into contiguous group. Cols belonging to
  // thread tid are
  //     thr_cols_seps[tid]:thr_col_seps[tid+1]-1
  // and similarly for rows.
  SetThreadIndexSeps(thr_row_seps, nthreads, m);
  SetThreadIndexSeps(thr_col_seps, nthreads, n);

#if 0
#ifndef NDEBUG
  printf("SetupParallel nnzs: ");
  for (size_t it = 0; it < thr_block_seps.size() - 1; it++) {
    size_t nnz = 0;
    for (size_t i = thr_block_seps[it]; i < thr_block_seps[it+1]; i++) {
      const Block<Real>& b = blocks[i];
      nnz += b.B ? b.m*b.n : b.rank*(b.m + b.n);
    }
    printf(" %ld", nnz);
  }
  printf("\n");
#endif
#endif
}
#endif

template<typename Real> void TypedHmat<Real>::ReorganizeMemory () {
  // Reallocate data into one large array to get contiguous memory.
  Real* old_blocks_mem = blocks_mem;
  try {
    // Since I'm allocating potentially a huge chunk, I want to safeguard it.
    blocks_mem = new Real[nnz];
  } catch (const bad_alloc&) {
    // Silently fail. We just won't reorganize memory. No big deal.
    blocks_mem = old_blocks_mem;
    return;
  }
  Real* mem = blocks_mem;
  for (size_t i = 0; i < blocks.size(); i++)
    mem += blocks[i].UseThisMemory(mem);
  if (old_blocks_mem) delete[] old_blocks_mem;
}

// xqi = x(q,:)
template<typename Real> template<typename T> inline void TypedHmat<Real>::
ApplyQ (const T* x, T* xqi, UInt ncol) const {
  // General comment about OpenMP blocks and function calls: If these MVP
  // routines are called from OpenMP threads -- i.e, from within a caller's
  // parallel block -- then the OpenMP stuff here is either nested if available
  // or simply ignored if not. Nonetheless, I've separated out the case nthreads
  // == 1 with the idea that I can avoid some needless OpenMP overhead. Not sure
  // if it's actually helpful.
  if (nthreads > 1) {
    omp_set_num_threads(nthreads);
    for (UInt ic = 0; ic < ncol; ic++) {
      UInt cos = ic*n;
#pragma omp parallel
      { UInt tid = omp_get_thread_num();
        for (UInt i = thr_col_seps[tid], N = thr_col_seps[tid+1];
             i < N; i++)
          xqi[cos + i] = (Real) x[cos + q[i]]; }
    }
  } else {
    for (UInt ic = 0; ic < ncol; ic++) {
      UInt cos = ic*n;
      for (UInt i = 0; i < n; i++)
        xqi[cos + i] = (Real) x[cos + q[i]];
    }
  }
}
  
// y(p,:) = yp
template<typename Real> template<typename T> inline void TypedHmat<Real>::
ApplyPt (const T* yp, T* y, UInt ncol) const {
  if (nthreads > 1) {
    omp_set_num_threads(nthreads);
#pragma omp parallel
    { UInt tid = omp_get_thread_num();
      Blint i0 = thr_row_seps[tid], i1 = thr_row_seps[tid+1];
      for (UInt ic = 0; ic < ncol; ic++) {
        UInt cos = ic*m;
        for (Blint j = i0; j < i1; j++)
          y[cos + p[j]] = (T) yp[cos + j];
      }}
  } else {
    for (UInt ic = 0; ic < ncol; ic++) {
      UInt cos = ic*m;
      for (Blint j = 0; j < m; j++)
        y[cos + p[j]] = (T) yp[cos + j];
    }
  }
}

// --> No (rs,cs)

template<typename Real> template<typename T> inline void TypedHmat<Real>::
PreMult (const T* x, UInt nthreads, UInt ncol) {
  //todo Possibly move col loop into parallel
  if (nthreads > 1) {
    if (do_perms) {
      // xq = x(q)
      for (UInt ic = 0; ic < ncol; ic++) {
        UInt cos = ic*n;
#pragma omp parallel
        { UInt tid = omp_get_thread_num();
          for (UInt i = thr_col_seps[tid], N = thr_col_seps[tid+1];
               i < N; i++)
            xq[cos + i] = (Real) x[cos + q[i]]; }
      }
    } else {
      // xq = x
      for (UInt ic = 0; ic < ncol; ic++) {
        UInt cos = ic*n;
#pragma omp parallel
        { UInt tid = omp_get_thread_num();
          for (UInt i = thr_col_seps[tid], N = thr_col_seps[tid+1];
               i < N; i++)
            xq[cos + i] = (Real) x[cos + i]; }
      }
    }
#pragma omp parallel
    { UInt tid = omp_get_thread_num();
      memset(ys + tid*ncol*m, 0, ncol*m*sizeof(Real)); }
  } else {
    if (do_perms) {
      // xq = x(q)
      for (UInt ic = 0; ic < ncol; ic++) {
        UInt cos = ic*n;
        for (UInt i = 0; i < n; i++)
          xq[cos + i] = (Real) x[cos + q[i]];
      }
    } else {
      // xq = x
      for (UInt i = 0, N = ncol*n; i < N; i++)
        xq[i] = (Real) x[i];
    }

    memset(ys, 0, ncol*m*sizeof(Real));
  }
}

template<typename Real> template<typename T> inline void TypedHmat<Real>::
PostMult (UInt nthreads, UInt ncol, T* y) {
  if (nthreads > 1) {
#pragma omp parallel
    { UInt tid = omp_get_thread_num();
      UInt i0 = thr_row_seps[tid], i1 = thr_row_seps[tid+1];
      memset(y + ncol*i0, 0, ncol*(i1 - i0)*sizeof(T)); }
#pragma omp parallel
    { UInt tid = omp_get_thread_num();
      Blint i0 = thr_row_seps[tid], i1 = thr_row_seps[tid+1];
      for (UInt ic = 0; ic < ncol; ic++) {
        UInt cos = ic*m;
        // Accumulate distributed y(:,ic)
        for (UInt i = 0; i < nthreads; i++) {
          int os = i*ncol*m;
          for (Blint j = i0; j < i1; j++) {
            // Pretty sure that compiler optimization should get rid of the
            // branch
            y[cos + (do_perms ? p[j] : j)] += (T) ys[cos + os + j];
          }
        }
      }}
  } else {
    memset(y, 0, ncol*m*sizeof(T));
    for (UInt ic = 0; ic < ncol; ic++) {
      UInt cos = ic*m;
      for (Blint j = 0; j < m; j++)
        y[cos + (do_perms ? p[j] : j)] += (T) ys[cos + j];
    }
  }
}

template<typename Real> template<typename T> void TypedHmat<Real>::
Mvp (const T* x, T* y, UInt ncol) throw (Exception) {
#ifdef TV
  struct timeval t1, t2, t3, t4;
  gettimeofday(&t1, 0);
#endif

  TypedHmat& h = *this;
  const bool transp = h.transp_mode;

  if (ncol > h.ncol) throw Exception("ncol > this->ncol");

  // If other OpenMP-parallelized routines have messed with this, we need to
  // reestablish what we want.
  if (nthreads > 1) omp_set_num_threads(nthreads);

  PreMult(x, nthreads, ncol);

#ifdef TV
  gettimeofday(&t3, 0);
#endif
  if (nthreads > 1) {
    if (ncol > 1) {
      // Only this is new
#pragma omp parallel
      { UInt tid = omp_get_thread_num();
        Real* yp = h.ys + tid*ncol*h.m;
        Real* wp = h.work + h.K*ncol*tid;
        for (Blint i = thr_block_seps[tid], N = thr_block_seps[tid+1];
             i < N; i++) {
          // ys += B*xq
          h.blocks[i].Mvp(transp, h.xq, yp, wp, h.m, h.n, ncol);
        }}
    } else {
#pragma omp parallel
      { UInt tid = omp_get_thread_num();
        Real* yp = h.ys + tid*h.m;
        Real* wp = h.work + h.K*tid;
        for (Blint i = thr_block_seps[tid], N = thr_block_seps[tid+1];
             i < N; i++) {
          // ys += B*xq
          h.blocks[i].Mvp(transp, h.xq, yp, wp);
        }}
    }
  } else {
    if (ncol > 1)
      for (UInt i = 0; i < h.nb; i++)
        h.blocks[i].Mvp(transp, h.xq, h.ys, h.work, h.m, h.n, ncol);
    else
      for (UInt i = 0; i < h.nb; i++)
        h.blocks[i].Mvp(transp, h.xq, h.ys, h.work);
  }
#ifdef TV
  gettimeofday(&t4, 0);
#endif
  
  PostMult(nthreads, ncol, y);

#ifdef TV
  gettimeofday(&t2, 0);
  printf("%f\n", difftime(t3, t4)/difftime(t1, t2));
#endif
}

// Swap a bunch of member data to implement the transpose MVP.
template<typename Real>
inline void TypedHmat<Real>::SwapForTranspose () {
  swap(m, n);
  swap(xq, ys);
  swap(p, q);
  swap(ip, iq);
  swap(rs, cs);
  swap(rs_idxs, cs_idxs);
  swap(irs, ics);
  swap(thr_row_seps, thr_col_seps);
  transp_mode = !transp_mode;
}

template<typename Real> template<typename T>
void TypedHmat<Real>::MvpT (const T* x, T* y, UInt ncol) throw (Exception) {
  SwapForTranspose();
  Mvp(x, y, ncol);
  SwapForTranspose();
}

// --> rs but no cs

// Inverse permutations of p and q
template<typename Real> inline void TypedHmat<Real>::InvPerms () {
  if (!ip && do_perms) {
    ip = new Blint[m];
    for (UInt i = 0; i < m; i++) ip[p[i]] = i;
    iq = new Blint[n];
    for (UInt i = 0; i < n; i++) iq[q[i]] = i;
  }
}

// I name and discuss variables in this routine in terms of rows, but it works
// on cols too.
inline void MapIdxs (
  const vector<Blint>& vrs, const Blint* ip, Blint*& rs, Blint m,
  Blint*& rs_idxs, bool do_e2b_map, bool do_perms, Blint* irs = 0) {
  // Lazily create array.
  if (!rs) rs = new Blint[m];
  if (!rs_idxs) rs_idxs = new Blint[m + 1];

  // vrs is the unordered set of the caller's row subset. Map vrs, which
  // applies to PBQ, to rs, which applies to B, using ip. (vcs -> cs using
  // iq.)
  UInt nrs = vrs.size();
  if (do_perms)
    for (UInt i = 0; i < nrs; i++) rs[i] = ip[vrs[i]];
  else
    memcpy(rs, &vrs[0], nrs*sizeof(Blint));

  // Inverse perm of rs if desired.
  if (irs) for (UInt i = 0; i < nrs; i++) irs[rs[i]] = i;

  // Sort rs for indexing in the next step.
  std::sort(rs, rs + nrs);

  if (!do_e2b_map) {
    // Index into rs. If i and k index the rows, say, then
    // rs(rs_idxs[i]:rs_idxs[k]-1) is the desired subset of i:k.
    Blint rp = -1;
    for (UInt i = 0; i < nrs; i++) {
      Blint r = rs[i];
      for (Blint j = rp + 1; j <= r; j++) rs_idxs[j] = i;
      rp = r;
    }
    // Set to rs.size() for the last part.
    for (Blint j = rp + 1; j <= m; j++) rs_idxs[j] = nrs;
  }
}

// Make use of rs
template<typename Real> template<typename T> inline void TypedHmat<Real>::
PostMult (UInt nthreads, UInt ncol, T* y, const vector<Blint>& vrs) {
  UInt nrs = vrs.size();
  for (UInt ic = 0; ic < ncol; ic++) {
    UInt cos = ic*m;
    if (nthreads == 1) {
      for (UInt j = 0; j < nrs; j++) {
        Blint r = vrs[j];
        y[cos + r] = (T) ys[cos + (do_perms ? ip[r] : r)];
      }
    } else {
      // Zero.
      for (UInt j = 0; j < nrs; j++)
        y[cos + vrs[j]] = 0.0;
      // Accumulate distributed y(:,ic).
      for (UInt i = 0; i < nthreads; i++) {
        int os = i*ncol*m;
        for (UInt j = 0; j < nrs; j++) {
          Blint r = vrs[j];
          y[cos + r] += (T) ys[cos + os + (do_perms ? ip[r] : r)];
        }
      }
    }
  }
}

template<typename Real> inline bool TypedHmat<Real>::
RIdxLimits (const Block<Real>& b, UInt nrs, Blint& rlo, Blint& rhi,
            Blint& nbrs) {
  rlo = rs_idxs[b.r0];
  if (rlo == (Blint) nrs) return true;
  rhi = rs_idxs[b.r0 + b.m];
  if (rhi == 0) return true;
  nbrs = rhi - rlo;
  return false;
}

template<typename Real> inline void TypedHmat<Real>::
RIdxLimitsNoTest (const Block<Real>& b, UInt nrs, Blint& rlo, Blint& rhi,
                  Blint& nbrs) {
  rlo = rs_idxs[b.r0];
  rhi = rs_idxs[b.r0 + b.m];
  nbrs = rhi - rlo;
}

void inline PushAndDoCapacity (vector<Blint>& v, UInt i) {
  UInt n = v.size();
  if (v.capacity() == n) v.reserve(2*n);
  v.push_back(i);
}

void inline Compact (vector<Blint>& v) {
  // Apparently, this is the standard method of compacting memory.
  vector<Blint>(v).swap(v);
}

template<typename Real>
inline void E2BMap::Init (const TypedHmat<Real>& hm) {
  if (is[0].size() == 0) {
    UInt m = hm.GetM(), n = hm.GetN();
    nb = hm.GetNbrBlocks();
    VecVecBlint& rs = is[0];
    VecVecBlint& cs = is[1];
    rs.resize(m);
    cs.resize(n);

    const vector<Block<Real> >& blocks = hm.GetBlocks();
    for (UInt i = 0; i < nb; i++) {
      const Block<Real>& b = blocks[i];
      for (Blint ic = 0; ic < b.n; ic++)
        PushAndDoCapacity(cs[b.c0 + ic], i);
      for (Blint ir = 0; ir < b.m; ir++)
        PushAndDoCapacity(rs[b.r0 + ir], i);
    }

    // Shrink memory to exactly what is needed.
    for (UInt i = 0; i < n; i++) Compact(cs[i]);
    for (UInt i = 0; i < m; i++) Compact(rs[i]);

    flag = 0;
    use_block = new FlagType[nb];
    memset(use_block, 0, nb*sizeof(FlagType));
  }
}

inline void E2BMap::InitFlag () {
  if (flag == (FlagType) -1) {
    flag = 0;
    memset(use_block, 0, nb*sizeof(FlagType));
  }
}

inline UInt SetIdxMapsE2B_GetIdx (const Blint* rs, UInt nrs, Blint r0) {
  UInt j;
  j = upper_bound(rs, rs + nrs, r0) - rs;
  if (j > 0 && rs[j-1] == r0) j--;
  return j;
}

template<typename Real>
inline void TypedHmat<Real>::SetIdxMapsE2B (UInt nrs, UInt ncs) {
  bool do_rs = nrs > 0, do_cs = ncs > 0;
  for (UInt i = 0; i < kblks; i++) {
    const Block<Real>& b = blocks[i];
    if (do_rs) {
      rs_idxs[b.r0] = SetIdxMapsE2B_GetIdx(rs, nrs, b.r0);
      rs_idxs[b.r0 + b.m] = SetIdxMapsE2B_GetIdx(rs, nrs, b.r0 + b.m);
    }
    if (do_cs) {
      cs_idxs[b.c0] = SetIdxMapsE2B_GetIdx(cs, ncs, b.c0);
      cs_idxs[b.c0 + b.n] = SetIdxMapsE2B_GetIdx(cs, ncs, b.c0 + b.n);
    }
  }
}

// Phrased in terms of columns, but equally valid for rows.
inline UInt E2BMap::
GetUsedBlocks (UInt cdim, const Blint* cs, UInt ncs, Blint* blks) {
  InitFlag();
  flag++;
  Blint kblks = 0;
  for (UInt ic = 0; ic < ncs; ic++) {
    const vector<Blint>& cv = is[cdim][cs[ic]];
    for (UInt i = 0, n = cv.size(); i < n; i++) {
      Blint cvi = cv[i];
      if (use_block[cvi] != flag) {
        use_block[cvi] = flag;
        blks[kblks] = cvi;
        kblks++;
      }
    }
  }

  return kblks;
}

template<typename Real> inline void TypedHmat<Real>::RGetUsedBlocks (UInt nrs) {
  if (!state_saved) {
    if (!blks) blks = new Blint[nb];

    if (do_e2b_map) {
      e2bm.Init(*this);
      kblks = e2bm.GetUsedBlocks(0, rs, nrs, blks);
      SetIdxMapsE2B(nrs, 0);
    } else {
      Blint rlo, rhi, nbrs;
      kblks = 0;
      for (UInt i = 0; i < nb; i++) {
        if (RIdxLimits(blocks[i], nrs, rlo, rhi, nbrs)) continue;
        if (nbrs > 0) {
          blks[kblks] = i;
          kblks++;
        }}
    }

    if (save_state) state_saved = true;

#if SORT_BLOCKS
    sort(blks, blks + kblks);
#endif
  }
}

template<typename Real> inline void TypedHmat<Real>::
RMvpBlock (UInt ncol, const Block<Real>& b, UInt nrs) {
  TypedHmat& h = *this;
  Blint rlo, rhi, nbrs;
  RIdxLimitsNoTest(b, nrs, rlo, rhi, nbrs);
  if (ncol == 1) {
    if (nbrs == b.m) b.Mvp(transp_mode, h.xq, h.ys, h.work);
    else b.RMvp(h.xq, h.ys, h.work, h.rs + rlo, nbrs);
  } else {
    for (UInt ic = 0; ic < ncol; ic++) {
      if (nbrs == b.m)
        b.Mvp(transp_mode, h.xq + ic*h.n, h.ys + ic*h.m, h.work);
      else
        b.RMvp(h.xq + ic*h.n, h.ys + ic*h.m, h.work,
               h.rs + rlo, nbrs);
    }
  }
}

template<typename Real> template<typename T> void TypedHmat<Real>::
RMvp (const T* x, T* y, UInt ncol, const vector<Blint>& vrs)
  throw (Exception)
{
  TypedHmat& h = *this;
  if (ncol > h.ncol) throw Exception("ncol > this->ncol");

  if (nthreads > 1) omp_set_num_threads(nthreads);
  InvPerms();
  if (!state_saved)
    MapIdxs(vrs, h.ip, h.rs, h.m, h.rs_idxs, h.do_e2b_map, h.do_perms);

  PreMult(x, nthreads, ncol);

  UInt nrs = vrs.size();
  Blint rlo, rhi, nbrs;
  if (nthreads == 1) {
    if (save_state || do_e2b_map) {
      RGetUsedBlocks(nrs);
      for (UInt i = 0; i < kblks; i++) {
        const Block<Real>& b = h.blocks[blks[i]];
        RMvpBlock(ncol, b, nrs);
      }
    } else {
      for (UInt i = 0; i < h.nb; i++) {
        const Block<Real>& b = h.blocks[i];
        if (RIdxLimits(b, nrs, rlo, rhi, nbrs)) continue;
        if (nbrs > 0) RMvpBlock(ncol, b, nrs);
      }
    }
  } else {
    RGetUsedBlocks(nrs);
#pragma omp parallel for private(rlo,rhi,nbrs) schedule(dynamic)
    for (UInt i = 0; i < kblks; i++) {
      const Block<Real>& b = h.blocks[blks[i]];
      RIdxLimitsNoTest(b, nrs, rlo, rhi, nbrs);
      UInt tid = omp_get_thread_num();
      if (ncol == 1) {
        if (nbrs == b.m)
          b.Mvp(transp_mode, h.xq, h.ys + tid*h.m, h.work + h.K*tid);
        else
          b.RMvp(h.xq, h.ys + tid*h.m, h.work + h.K*tid, h.rs + rlo, nbrs);
      }
      else {
        if (nbrs == b.m)
          for (UInt ic = 0; ic < ncol; ic++)
            b.Mvp(transp_mode, h.xq + ic*h.n,
                  h.ys + (ic + tid*ncol)*h.m,
                  h.work + h.K*tid);
        else
          for (UInt ic = 0; ic < ncol; ic++) {
            b.RMvp(h.xq + ic*h.n,
                   h.ys + (ic + tid*ncol)*h.m,
                   h.work + h.K*tid,
                   h.rs + rlo, nbrs);
          }}}
  }
  
  PostMult(nthreads, ncol, y, vrs);
}

// --> cs but no rs

template<typename Real> inline bool TypedHmat<Real>::
CIdxLimits (const Block<Real>& b, UInt ncs, Blint& clo, Blint& chi,
            Blint& nbcs) {
  clo = cs_idxs[b.c0];
  if (clo == (Blint) ncs) return true;
  chi = cs_idxs[b.c0 + b.n];
  if (chi == 0) return true;
  nbcs = chi - clo;
  return false;
}

template<typename Real> inline void TypedHmat<Real>::
CIdxLimitsNoTest (const Block<Real>& b, UInt crs, Blint& clo, Blint& chi,
                  Blint& nbcs) {
  clo = cs_idxs[b.c0];
  chi = cs_idxs[b.c0 + b.n];
  nbcs = chi - clo;
}

template<typename Real> inline void TypedHmat<Real>::CGetUsedBlocks (UInt ncs) {
  if (!state_saved) {
    if (!blks) blks = new Blint[nb];

    if (do_e2b_map) {
      e2bm.Init(*this);
      kblks = e2bm.GetUsedBlocks(1, cs, ncs, blks);
      SetIdxMapsE2B(0, ncs);
    } else {
      Blint clo, chi, nbcs;
      kblks = 0;
      for (UInt i = 0; i < nb; i++) {
        if (CIdxLimits(blocks[i], ncs, clo, chi, nbcs)) continue;
        if (nbcs > 0) {
          blks[kblks] = i;
          kblks++;
        }}
    }

    if (save_state) state_saved = true;
#if SORT_BLOCKS
    sort(blks, blks + kblks);
#endif
  }
}

// Make use of cs
template<typename Real> template<typename T> inline void TypedHmat<Real>::
PreMult (const T* x, UInt nthreads, UInt ncol, const vector<Blint>& vcs) {
  UInt ncs = vcs.size();
  // xq(iq(cs)) = x(cs)
  for (UInt ic = 0; ic < ncol; ic++) {
    UInt cos = ic*n;
    for (UInt i = 0; i < ncs; i++) {
      Blint c = vcs[i];
      xq[cos + (do_perms ? iq[c] : c)] = (Real) x[cos + c];
    }
  }
  memset(ys, 0, nthreads*ncol*m*sizeof(Real));
}

template<typename Real> inline void TypedHmat<Real>::
CMvpBlock (UInt ncol, const Block<Real>& b, UInt ncs) {
  TypedHmat& h = *this;
  Blint clo, chi, nbcs;
  CIdxLimitsNoTest(b, ncs, clo, chi, nbcs);
  if (ncol == 1) {
    if (nbcs == b.n) b.Mvp(transp_mode, h.xq, h.ys, h.work);
    else b.CMvp(h.xq, h.ys, h.work, h.cs + clo, nbcs);
  } else {
    for (UInt ic = 0; ic < ncol; ic++) {
      if (nbcs == b.n)
        b.Mvp(transp_mode, h.xq + ic*h.n, h.ys + ic*h.m, h.work);
      else
        b.CMvp(h.xq + ic*h.n, h.ys + ic*h.m, h.work,
               h.cs + clo, nbcs);
    }
  }
}

template<typename Real> template<typename T> void TypedHmat<Real>::
CMvp (const T* x, T* y, UInt ncol, const vector<Blint>& vcs)
  throw (Exception)
{
  TypedHmat& h = *this;
  if (ncol > h.ncol) throw Exception("ncol > this->ncol");

  if (nthreads > 1) omp_set_num_threads(nthreads);
  InvPerms();
  if (!state_saved)
    MapIdxs(vcs, h.iq, h.cs, h.n, h.cs_idxs, h.do_e2b_map, h.do_perms);

  PreMult(x, nthreads, ncol, vcs);

  UInt ncs = vcs.size();
  Blint clo, chi, nbcs;
  if (nthreads == 1) {
    if (save_state || do_e2b_map) {
      CGetUsedBlocks(ncs);
      for (UInt i = 0; i < kblks; i++) {
        const Block<Real>& b = h.blocks[blks[i]];
        CMvpBlock(ncol, b, ncs);
      }
    } else {
      for (UInt i = 0; i < h.nb; i++) {
        const Block<Real>& b = h.blocks[i];
        if (CIdxLimits(b, ncs, clo, chi, nbcs)) continue;
        if (nbcs > 0) CMvpBlock(ncol, b, ncs);
      }
    }
  } else {
    CGetUsedBlocks(ncs);
#pragma omp parallel for private(clo,chi,nbcs) schedule(dynamic)
    for (UInt i = 0; i < kblks; i++) {
      const Block<Real>& b = h.blocks[blks[i]];
      CIdxLimitsNoTest(b, ncs, clo, chi, nbcs);
      UInt tid = omp_get_thread_num();
      if (ncol == 1) {
        if (nbcs == b.n)
          b.Mvp(transp_mode, h.xq, h.ys + tid*h.m, h.work + h.K*tid);
        else
          b.CMvp(h.xq, h.ys + tid*h.m, h.work + h.K*tid, h.cs + clo, nbcs);
      }
      else {
        if (nbcs == b.n)
          for (UInt ic = 0; ic < ncol; ic++)
            b.Mvp(transp_mode, h.xq + ic*h.n,
                  h.ys + (ic + tid*ncol)*h.m,
                  h.work + h.K*tid);
        else
          for (UInt ic = 0; ic < ncol; ic++) {
            b.CMvp(h.xq + ic*h.n,
                   h.ys + (ic + tid*ncol)*h.m,
                   h.work + h.K*tid,
                   h.cs + clo, nbcs);
          }}}
  }
  
  PostMult(nthreads, ncol, y);
}

// --> (rs,cs)

// Make use of rs and cs
template<typename Real> template<typename T> inline void TypedHmat<Real>::
PreMult (const T* x, UInt nthreads, UInt ncol, const vector<Blint>& vrs,
         const vector<Blint>& vcs) {
  UInt ncs = vcs.size();
  // xq(iq(cs)) = x(cs)
  for (UInt ic = 0; ic < ncol; ic++) {
    UInt cos = ic*n;
    for (UInt i = 0; i < ncs; i++) {
      Blint c = vcs[i];
      xq[cos + (do_perms ? iq[c] : c)] = (Real) x[cos + c];
    }
  }
  //memset(ys, 0, nthreads*ncol*m*sizeof(Real));
  UInt os = 0;
  for (UInt i = 0; i < nthreads; i++) {
    for (UInt ic = 0; ic < ncol; ic++) {
      for (UInt ir = 0, nrs = vrs.size(); ir < nrs; ir++) {
        Blint r = vrs[ir];
        ys[os + (do_perms ? ip[r] : r)] = 0.0;
      }
      os += m;
    }
  }
}

template<typename Real> inline bool TypedHmat<Real>::
IdxLimits (const Block<Real>& b, UInt nrs, UInt ncs, Blint& rlo, Blint& rhi,
           Blint& clo, Blint& chi, Blint& nbrs, Blint& nbcs) {
  rlo = rs_idxs[b.r0];
  if (rlo == (Blint) nrs) return true;
  rhi = rs_idxs[b.r0 + b.m];
  if (rhi == 0) return true;
  clo = cs_idxs[b.c0];
  if (clo == (Blint) ncs) return true;
  chi = cs_idxs[b.c0 + b.n];
  if (chi == 0) return true;
  nbrs = rhi - rlo;
  nbcs = chi - clo;
  return false;
}

template<typename Real> inline void TypedHmat<Real>::
IdxLimitsNoTest (const Block<Real>& b, UInt nrs, UInt ncs, Blint& rlo,
                 Blint& rhi, Blint& clo, Blint& chi, Blint& nbrs, Blint& nbcs) {
  rlo = rs_idxs[b.r0];
  rhi = rs_idxs[b.r0 + b.m];
  clo = cs_idxs[b.c0];
  chi = cs_idxs[b.c0 + b.n];
  nbrs = rhi - rlo;
  nbcs = chi - clo;
}

template<typename Real> inline void TypedHmat<Real>::
GetUsedBlocks (UInt nrs, UInt ncs) {
  if (!state_saved) {
    if (!blks) blks = new Blint[nb];
    
    if (do_e2b_map) {
      e2bm.Init(*this);
#ifdef TV
      struct timeval t1, t2, t3;
      gettimeofday(&t1, 0);
#endif
      kblks = GetUsedBlocksE2B(nrs, ncs);
#ifdef TV
      gettimeofday(&t2, 0);
#endif
      SetIdxMapsE2B(nrs, ncs);
#ifdef TV
      gettimeofday(&t3, 0);
      printf("GUB: %f\n", difftime(t1, t2)/difftime(t1, t3));
#endif
    } else {
      // Find out which blocks are used
      Blint rlo, rhi, clo, chi, nbrs, nbcs;
      kblks = 0;
      for (UInt i = 0; i < nb; i++) {
        if (IdxLimits(blocks[i], nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs))
          continue;
        if (nbrs > 0 && nbcs > 0) {
          blks[kblks] = i;
          kblks++;
        }
      }
    }
#if SORT_BLOCKS
    sort(blks, blks + kblks);
#endif
    if (save_state) state_saved = true;
  }
}

UInt E2BMap::GetUsedBlocks (UInt cdim, const Blint* cs, UInt ncs,
                            UInt rdim, const Blint* rs, UInt nrs, Blint* blks) {
  // Get all blocks needed for the column subset.
  InitFlag();
  flag++;
  for (UInt ic = 0; ic < ncs; ic++) {
    const vector<Blint>& cv = is[cdim][cs[ic]];
    for (UInt i = 0, n = cv.size(); i < n; i++)
      use_block[cv[i]] = flag;     // This block might be needed.
  }
  // Intersect with those needed for the row subset.
  Blint kblks = 0;
  for (UInt ir = 0; ir < nrs; ir++) {
    const vector<Blint>& rv = is[rdim][rs[ir]];
    for (UInt i = 0, n = rv.size(); i < n; i++) {
      UInt rvi = rv[i];
      if (use_block[rvi] == flag) {
        use_block[rvi] = flag - 1; // Insert block only once.
        blks[kblks] = rvi;         // Insert block: it's definitely needed.
        kblks++;
      }
    }
  }

  return kblks;
}                                                         

template<typename Real> inline UInt TypedHmat<Real>::
GetUsedBlocksE2B (UInt nrs, UInt ncs) {
  if (ncs <= nrs)
    return e2bm.GetUsedBlocks(1, cs, ncs, 0, rs, nrs, blks);
  else
    return e2bm.GetUsedBlocks(0, rs, nrs, 1, cs, ncs, blks);
}

template<typename Real> inline void TypedHmat<Real>::
MvpBlock (UInt ncol, const Block<Real>& b, UInt nrs, UInt ncs) {
  TypedHmat& h = *this;
  Blint rlo, rhi, nbrs, clo, chi, nbcs;
  IdxLimitsNoTest(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs);
  if (ncol == 1) {
    if (nbrs == b.m && nbcs == b.n) b.Mvp(transp_mode, h.xq, h.ys, h.work);
    else b.Mvp(h.xq, h.ys, h.work, h.rs + rlo, nbrs, h.cs + clo, nbcs);
  } else {
    for (UInt ic = 0; ic < ncol; ic++) {
      if (nbrs == b.m && nbcs == b.n)
        b.Mvp(transp_mode, h.xq + ic*h.n, h.ys + ic*h.m, h.work);
      else
        b.Mvp(h.xq + ic*h.n, h.ys + ic*h.m, h.work,
              h.rs + rlo, nbrs, h.cs + clo, nbcs);
    }
  }
}

template<typename Real> template<typename T> void TypedHmat<Real>::
Mvp (const T* x, T* y, UInt ncol, const vector<Blint>& vrs,
     const vector<Blint>& vcs)
  throw (Exception)
{
#ifdef TV
  struct timeval t1, t2, t3, t4, t5, t6;
  double et = 0.0;
  gettimeofday(&t1, 0);
#endif

  TypedHmat& h = *this;
  if (ncol > h.ncol) throw Exception("ncol > this->ncol");

  if (nthreads > 1) omp_set_num_threads(nthreads);
  InvPerms();
  if (!state_saved) {
    MapIdxs(vrs, h.ip, h.rs, h.m, rs_idxs, h.do_e2b_map, h.do_perms);
    MapIdxs(vcs, h.iq, h.cs, h.n, cs_idxs, h.do_e2b_map, h.do_perms);
  }

  PreMult(x, nthreads, ncol, vrs, vcs);

  UInt nrs = vrs.size(), ncs = vcs.size();
  Blint rlo, rhi, nbrs, clo, chi, nbcs;
#ifdef TV
  gettimeofday(&t3, 0);
#endif
  if (nthreads == 1) {
    if (save_state || do_e2b_map) {
#ifdef TV
      gettimeofday(&t5, 0);
#endif
      GetUsedBlocks(nrs, ncs);
#ifdef TV
      gettimeofday(&t6, 0);
#endif
      for (UInt i = 0; i < kblks; i++) {
        const Block<Real>& b = h.blocks[blks[i]];
        MvpBlock(ncol, b, nrs, ncs);
      }
    } else {
      for (UInt i = 0; i < h.nb; i++) {
        const Block<Real>& b = h.blocks[i];
        if (IdxLimits(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs)) continue;
        if (nbrs > 0 && nbcs > 0) MvpBlock(ncol, b, nrs, ncs);
      }
    }
  } else {
    GetUsedBlocks(nrs, ncs);
#pragma omp parallel for private(rlo,rhi,clo,chi,nbrs,nbcs) schedule(dynamic)
    for (UInt i = 0; i < kblks; i++) {
      const Block<Real>& b = h.blocks[blks[i]];
      IdxLimitsNoTest(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs);
      UInt tid = omp_get_thread_num();
      if (ncol == 1) {
        if (nbrs == b.m && nbcs == b.n)
          b.Mvp(transp_mode, h.xq, h.ys + tid*h.m, h.work + h.K*tid);
        else
          b.Mvp(h.xq, h.ys + tid*h.m, h.work + h.K*tid,
                h.rs + rlo, nbrs, h.cs + clo, nbcs);
      } else {
        for (UInt ic = 0; ic < ncol; ic++) {
          b.Mvp(h.xq + ic*h.n,
                h.ys + (ic + tid*ncol)*h.m,
                h.work + h.K*tid,
                h.rs + rlo, nbrs, h.cs + clo, nbcs);
        }}}
  }
#ifdef TV
  gettimeofday(&t4, 0);
#endif

  PostMult(nthreads, ncol, y, vrs);
#ifdef TV
  gettimeofday(&t2, 0);
  if (nthreads == 1 && (save_state || do_e2b_map))
    printf("%f %f\n", (difftime(t3, t4) - difftime(t5, t6))/difftime(t1, t2),
           difftime(t3, t4)/difftime(t1, t2));
  else
    printf("%f\n", difftime(t3, t4)/difftime(t1, t2));
#endif
}

// --> Other stuff

template<typename Real> UInt TypedHmat<Real>::
GetNnz (const vector<Blint>& vrs, const vector<Blint>& vcs) {
  TypedHmat& h = *this;
  InvPerms();
  if (!state_saved) {
    MapIdxs(vrs, h.ip, h.rs, h.m, h.rs_idxs, h.do_e2b_map, h.do_perms);
    MapIdxs(vcs, h.iq, h.cs, h.n, h.cs_idxs, h.do_e2b_map, h.do_perms);
  }

  UInt nnz = 0;
  UInt nrs = vrs.size(), ncs = vcs.size();  
  if (save_state) {
    GetUsedBlocks(nrs, ncs);
    for (UInt i = 0; i < kblks; i++) {
      const Block<Real>& b = h.blocks[blks[i]];
      Blint rlo, rhi, nbrs, clo, chi, nbcs;
      IdxLimitsNoTest(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs);
      nnz += b.Nnz(nbrs, nbcs);
    }
  } else {
    for (UInt i = 0; i < nb; i++) {
      const Block<Real>& b = h.blocks[i];
      Blint rlo, rhi, nbrs, clo, chi, nbcs;
      if (IdxLimits(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs)) continue;
      if (nbrs > 0 && nbcs > 0) nnz += b.Nnz(nbrs, nbcs);
    }
  }
  return nnz;
}

template<typename Real> double TypedHmat<Real>::NormFrobenius2 () const {
  const TypedHmat& h = *this;
  double nf2 = 0.0;
  if (nthreads > 1) {
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
    for (UInt i = 0; i < h.nb; i++) {
      double bnf2 = blocks[i].NormFrobenius2();
#pragma omp critical (cMP)
      { nf2 += bnf2; }
    }
  } else {
    for (UInt i = 0; i < h.nb; i++) nf2 += blocks[i].NormFrobenius2();
  }

  return nf2;
}

template<typename Real> double TypedHmat<Real>::NormOne () const {
  const TypedHmat& h = *this;
  if (nthreads > 1) omp_set_num_threads(nthreads);

  vector<double> col_norm1s(n, 0.0);
  for (UInt i = 0; i < h.nb; i++)
    blocks[i].AccumColNorm1s(work, col_norm1s);
  return *std::max_element(col_norm1s.begin(), col_norm1s.end());
}

template<typename Real> template<typename T> void TypedHmat<Real>::
Extract (const vector<Blint>& vrs, const vector<Blint>& vcs, T* es) {
  TypedHmat& h = *this;
  UInt nrs = vrs.size(), ncs = vcs.size();

  if (nthreads > 1) omp_set_num_threads(nthreads);
  InvPerms();
  if (!irs) {
    irs = new Blint[h.m];
    ics = new Blint[h.n];
  }
  if (!state_saved) {
    MapIdxs(vrs, h.ip, h.rs, h.m, rs_idxs, h.do_e2b_map, h.do_perms, irs);
    MapIdxs(vcs, h.iq, h.cs, h.n, cs_idxs, h.do_e2b_map, h.do_perms, ics);
  }

  Blint rlo, rhi, nbrs, clo, chi, nbcs;
  if (true || nthreads == 1) {
    if (save_state || do_e2b_map) {
      GetUsedBlocks(nrs, ncs);
      for (UInt i = 0; i < kblks; i++) {
        const Block<Real>& b = h.blocks[blks[i]];
        ExtractBlock(b, nrs, ncs, es);
      }
    } else {
      for (UInt i = 0; i < h.nb; i++) {
        const Block<Real>& b = h.blocks[i];
        if (IdxLimits(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs)) continue;
        if (nbrs > 0 && nbcs > 0) ExtractBlock(b, nrs, ncs, es);
      }
    }
  } else {
    // Not implemented yet. I'm guessing multithreading this routine won't be
    // useful.
  }
}

// In Hmat rather than Block because I want to access both global and local
// (rs,cs) information simultaneously.
template<typename Real> template<typename T> inline void TypedHmat<Real>::
ExtractBlock (const Block<Real>& b, Blint nrs, Blint ncs, T* es) {
  TypedHmat& h = *this;
  Blint rlo, rhi, nbrs, clo, chi, nbcs;
  IdxLimitsNoTest(b, nrs, ncs, rlo, rhi, clo, chi, nbrs, nbcs);
  if (b.B) {
    for (Blint j = 0; j < nbcs; j++) {
      Blint c = h.cs[clo + j];
      Blint cos = nrs*h.ics[c];
      const Real* B = b.B + (c - b.c0)*b.m;
      for (Blint i = 0; i < nbrs; i++) {
        Blint r = h.rs[rlo + i];
        es[cos + h.irs[r]] = (T) B[r - b.r0];
      }
    }
  } else {
    for (Blint j = 0; j < nbcs; j++) {
      Blint c = h.cs[clo + j];
      Blint cos = nrs*h.ics[c];
      const Real* Vt = b.Vt + b.rank*(c - b.c0);
      for (Blint i = 0; i < nbrs; i++) {
        Blint r = h.rs[rlo + i];
        // d = U(r,:)*Vt(:,c)
        T d = 0.0;
        const Real* U = b.U;
        for (Blint k = 0; k < b.rank; k++) {
          d += (T) U[r - b.r0]*Vt[k];
          U += b.m;
        }
        es[cos + h.irs[r]] = d;
      }}
  }
}

template<typename Real> template<typename T> void TypedHmat<Real>::
FullBlocksIJS (vector<Blint>& I, vector<Blint>& J, vector<T>& S, T cutoff)
  const
{
  const TypedHmat& h = *this;
  // Loop through twice to get cutoff and nnz for I,J,S. I don't want to use a
  // conservative estimate like nnz of this H-matrix, as that could be
  // disastrous. In general, this implementation doesn't care about iterating
  // through memory but rather using new memory.

  // Max magnitude.
  double max_mag = 0.0;
  if (cutoff > 0.0) {
    for (UInt i = 0; i < h.nb; i++)
      if (h.blocks[i].B) {
        double mm = h.blocks[i].MaxMag();
        if (mm > max_mag) max_mag = mm;
      }
  }
  cutoff *= max_mag;

  // Space to allocate.
  size_t idx = 0;
  for (UInt i = 0; i < h.nb; i++)
    if (h.blocks[i].B) {
      const Block<Real>& b = h.blocks[i];
      if (cutoff > 0.0) {
        for (int i = 0, N = b.m*b.n; i < N; i++)
          if (fabs(b.B[i]) >= cutoff) idx++;
      } else {
        idx += b.m*b.n;
      }
    }
  I.resize(idx); J.resize(idx); S.resize(idx);

  // Now grab each full-block block.
  idx = 0;
  for (UInt i = 0; i < h.nb; i++)
    if (h.blocks[i].B) {
      const Block<Real>& b = h.blocks[i];
      for (Blint k = 0, c = 0; k < b.n; k++)
        for (Blint j = 0; j < b.m; j++, c++) {
          if (cutoff > 0.0 && fabs(b.B[c]) < cutoff) continue;
          if (do_perms) {
            I[idx] = h.p[b.r0 + j] + 1;
            J[idx] = h.q[b.c0 + k] + 1;
          } else {
            I[idx] = b.r0 + j + 1;
            J[idx] = b.c0 + k + 1;
          }
          S[idx] = b.B[c];
          idx++;
        }
    }
}

#undef DELETE
}
