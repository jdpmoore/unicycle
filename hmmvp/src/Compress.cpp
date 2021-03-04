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

#include <stdio.h>
#include <math.h>
#include <sstream>
#include "util/include/Mpi.hpp"
#include "util/include/OpenMP.hpp"
#include "util/include/Util.hpp"
#include "util/include/CodeAnalysis.hpp"
#include "hmmvp/include/Hmat.hpp"
#include "hmmvp/include/HmatIo.hpp"
#include "Hd_pri.hpp"
#include "Compress_pri.hpp"

namespace hmmvp {
using namespace util;

// Pack and unpack basic objects into/from ByteBuffers.
template<typename T>
inline void PackVec (mpi::ByteBufferWriter& bw, const vector<T>& v) {
  bw.WriteScalar((size_t) v.size());
  bw.Write(&v[0], v.size());
}

template<typename T>
inline void UnpackVec (mpi::ByteBufferReader& br, vector<T>& v) {
  v.resize(br.ViewScalarAndAdvance<size_t>());
  br.Read(&v[0], v.size());
}

template<typename T>
inline void PackMatrix (mpi::ByteBufferWriter& bw, const Matrix<T>& A) {
  bw.WriteScalar((size_t) A.Size(1));
  bw.WriteScalar((size_t) A.Size(2));
  bw.Write(A.GetPtr(), A.Size());
}

template<typename T>
inline void UnpackMatrix (mpi::ByteBufferReader& br, Matrix<T>& A) {
  size_t m = br.ViewScalarAndAdvance<size_t>();
  size_t n = br.ViewScalarAndAdvance<size_t>();
  A.Resize(m, n);
  br.Read(A.GetPtr(), m*n);
}

const double LraOptions::qr_alpha = 0.5;

MatBlock::MatBlock () : r0(0), m(0), c0(0), n(0) {}

MatBlock::MatBlock (UInt r0_, UInt m_, UInt c0_, UInt n_)
  : r0(r0_), m(m_), c0(c0_), n(n_)
{
  assert(r0 >= 0 && c0 >= 0 && m >= 1 && n >= 1);
}

MatBlock& MatBlock::operator= (const hd::Block& b) {
  r0 = b.r0 - 1; m = b.m; c0 = b.c0 - 1; n = b.n;
  return *this;
}
  
template<typename T> void TypedLraBlock<T>::WriteToFile (FILE* fid) const {
  const MatBlock& mb = GetBlock();
  Blint n;
  n = mb.r0; write(&n, 1, fid);
  n = mb.m;  write(&n, 1, fid);
  n = mb.c0; write(&n, 1, fid);
  n = mb.n;  write(&n, 1, fid);
  if (HaveB()) {
    char s = 'B'; write(&s, 1, fid);
    write(B().GetPtr(), B().Size(), fid);
  } else {
    char s = 'U'; write(&s, 1, fid);
    n = U().Size(2); write(&n, 1, fid);
    write(U().GetPtr(), U().Size(), fid);
    Matrix<T> Vt;
    Transpose(V(), Vt);
    write(Vt.GetPtr(), Vt.Size(), fid);
  }
}

template<typename T> void TypedLraBlock<T>::
WriteInMemory (HmatConstructor& hmc, size_t bi) const {
  const MatBlock& mb = GetBlock();
  if (HaveB())
    hmc.SetBlockB(bi, mb.r0, mb.c0, B());
  else
    hmc.SetBlockUV(bi, mb.r0, mb.c0, U(), V());
}

template<typename T>
void TypedLraBlock<T>::Pack (mpi::ByteBufferWriter& bw) const {
  if (HaveB()) {
    bw.WriteScalar((char) 'B');
    PackMatrix(bw, B());
  } else {
    bw.WriteScalar((char) 'U');
    PackMatrix(bw, U());
    PackMatrix(bw, V());
  }
}

template<typename T> void TypedLraBlock<T>::Unpack (mpi::ByteBufferReader& br) {
  char code = br.ViewScalarAndAdvance<char>();
  if (code == 'B') {
    UnpackMatrix(br, B());
  } else {
    UnpackMatrix(br, U());
    UnpackMatrix(br, V());
  }
}

LraBlock*
NewLraBlock (const MatBlock& b, UInt realp, bool is_mrem, double tol) {
  if (realp == 1) return new TypedLraBlock<float>(b, is_mrem, tol);
  else return new TypedLraBlock<double>(b, is_mrem, tol);
}

// -----------------------------------------------------------------------------
// Low-rank approximation. No need to carry over all of the options in
// hm('Compress').

template<typename T>
class FullBMa : public MatrixAccessor {
public:
  FullBMa(const Matrix<T>* B) : _B(B) {}
  virtual bool Call(const CompressBlockInfo& cbi, const MatBlock& blk,
                    const vector<UInt>* rs, const vector<UInt>* cs, double* B)
    const;
private:
  const Matrix<T>* _B;
};

template<typename T>
bool FullBMa<T>::Call (const CompressBlockInfo& cbi, const MatBlock& blk,
                       const vector<UInt>* rs, const vector<UInt>* cs,
                       double* B) const {
  if (!rs && !cs) {
    Matrix<double> wB;
    wB.SetPtr(blk.m, blk.n, B);
    wB = *_B;
    return true;
  }
  if (!rs) {
    for (UInt ic = 0, nc = cs->size(), k = 0; ic < nc; ic++) {
      const T* p_B = _B->GetPtr() + ((*cs)[ic] - blk.c0) * blk.m;
      for (UInt ir = 0; ir < blk.m; ir++, k++) B[k] = (double) p_B[ir];
    }
    return true;
  }
  if (!cs) {
    const T* p_B = _B->GetPtr();
    UInt nr = rs->size();
    for (UInt ic = 0, k = 0; ic < blk.n; ic++) {
      for (UInt ir = 0; ir < nr; ir++, k++)
        B[k] = (double) p_B[(*rs)[ir] - blk.r0];
      p_B += blk.m;
    }
    return true;
  }
  // Should never have both rs and cs.
  assert(false);
  return false;
}

template<typename T>
void FillBFromGh (AcaGfHolder* gh, const CompressBlockInfo& cbi,
                  const MatBlock& blk, MatrixAccessor& ma, Matrix<T>& B) {
  B.Resize(blk.m, blk.n);
  T* pB = B.GetPtr();
  vector<double> data(blk.m);
  vector<UInt> col(1, blk.c0);
  for (size_t i = 0; i < blk.n; i++) {
    gh->Call(cbi, ma, NULL, &col, &data[0], false);
    for (size_t j = 0; j < blk.m; j++) pB[j] = (T) data[j];
    pB += blk.m;
    col[0]++;
  }
}

template<typename T, typename OT>
int ApproxByLowRank (const LraOptions& opts, TypedLraBlock<T>& b,
                     MatrixAccessor& ma, const TypedLraBlock<OT>* ob)
  throw (OutOfMemoryException, UserReqException)
{
  const MatBlock& blk = b.GetBlock();

  int ret = 0;
  if (ob) {
    // IsCompatible guarantees that we're both using the same of MREM and BREM.
    double old_tol = ob->GetTol();
    double new_tol = b.GetTol();
    // ||B_single - B_double||_F ~= eps(single) sqrt(m n)
    bool trust_ob = ob->GetPrec() == 2
      || new_tol >= numeric_limits<OT>::epsilon() *
      sqrt((double) blk.m * blk.n);
    if (trust_ob) {
      if (new_tol < old_tol && ob->HaveB()) {
        // We need a full B, and we can use the old one.
        b.B() = ob->B();
        return 1;
      } else if (new_tol > old_tol && !ob->HaveB()) {
        // We can use the old and U and V.
        b.U() = ob->U();
        b.V() = ob->V();
        if (b.NcolUV() > 1) {
          // We have ||C - A|| <= old_tol and want ||C - B|| <= new_tol. Hence
          //     ||C - B|| = ||(C - A) + (A - B)|| <= ||C - A|| + ||A - B||
          //       <= old_tol + tol <= new_tol
          // and so tol = new_tol - old_tol. (If BREM is being used, then just
          // multiply all terms by ||B||, and then drop it to keep tol in
          // block-wise relative form.)
          double tol = new_tol - old_tol;
          b.U() = ob->U();
          b.V() = ob->V();
          CompressQr(b.U(), b.V(), !opts.is_mrem, tol, opts.allow_0rank);
        }
        return 2;
      }
    } else {
      ob = NULL;
      ret = 100;
    }
  }

  if (blk.m == 1 || blk.n == 1 || (blk.m <= 2 && blk.n <= 2)) {
    // Too small, so just grab B.
    b.U().Resize(0);
    b.V().Resize(0);
    if (ob && ob->HaveB()) {
      ret = 10;
      b.B() = ob->B();
      return 1;
    } else {
      Matrix<double> B(blk.m, blk.n);
      ma.Call(CompressBlockInfo(blk.m, blk.n, opts.is_mrem, b.GetTol()), blk,
              NULL, NULL, B.GetPtr());
      b.B() = B;
      return ret + 0;
    }
  }

  LraOptions::LraMethod method = opts.method;
  if (std::min(blk.m, blk.n) < opts.min_aca_size) {
    // It's fine if the other dimension is very large, as we use the thin SVD.
    method = LraOptions::lra_svd;
  }

  switch (method) {
  case LraOptions::lra_svd: {
    Matrix<T> s, Vt;
    Matrix<double> B;
    if (ob && ob->HaveB()) {
      ret = 10;
      Svd(ob->B(), b.U(), s, Vt);
    } else {
      B.Resize(blk.m, blk.n);
      ma.Call(CompressBlockInfo(blk.m, blk.n, opts.is_mrem, b.GetTol()), blk,
              NULL, NULL, B.GetPtr());
      Svd(B, b.U(), s, Vt);
    }
    UInt rank = SvdChooseRank(s, !opts.is_mrem, b.GetTol(), opts.allow_0rank);
    if ((blk.m + blk.n) * rank < blk.m*blk.n) {
      T* pU = b.U().Reshape(blk.m, rank).GetPtr();
      for (UInt j = 0, k = 0; j < rank; j++) {
        T sj = s(j+1);
        for (UInt i = 0; i < blk.m; i++, k++) pU[k] = pU[k] * sj;
      }
      T* pV = b.V().Resize(blk.n, rank).GetPtr();
      for (UInt j = 0, k = 0; j < rank; j++)
        for (UInt i = 0; i < blk.n; i++, k++) pV[k] = Vt(j+1, i+1);
    } else {
      // Better just to keep the whole block.
      b.U().Resize(0);
      b.V().Resize(0);
      if (ob && ob->HaveB()) b.B() = ob->B();
      else b.B() = B;
    }
    return ret + 3;
  }

  case LraOptions::lra_aca: {
    MatrixAccessor* pma = &ma;
    bool ma_alloced = false;
    // Make use of old data, if available.
    if (ob) {
      if (ob->HaveB()) {
        ret = 20;
        ma_alloced = true;
        pma = new FullBMa<OT>(&ob->B());
      } else {
        ret = 30;
        b.U() = ob->U();
        b.V() = ob->V();
      }
    }

    // alpha splits up the error between Aca and CompressQr.
    double alpha = opts.call_CompressQr ? opts.qr_alpha : 1;
    // In addition, we set this first error to ~1/10 of the requested error
    // because ACA sometimes is not accurate enough by ~<= 10x. CompressQr
    // generally makes up for this. (aca_factor is nominally 1/10, but can have
    // another value.)
    double tol = alpha * opts.aca_tol_factor * b.GetTol();
    // -1.0 forces the first irow to be accepted (unless ob).
    double scale = -1.0;
    // Default ACA requests redundant Green's function evaluations. For matrices
    // whose approx rank is small compared to the matrix size, redundancy
    // doesn't matter, and if the GF eval is decently fast, it's much more
    // efficient than avoiding redundancy. Try to come up with a decent rule for
    // deciding when to avoid redundancy.
    bool avoid_redundant_gf_evals =
      // No reason to avoid calls if we're using a FullBMa.
      !ma_alloced
      && (opts.avoid_redundant_gf_evals || max(blk.m, blk.n) <= 2 << 10);
    AcaGfHolder* gh = NULL;
    if (avoid_redundant_gf_evals) gh = new AcaGfHolder(&blk);
    Aca(blk, *pma, scale, !opts.is_mrem, tol, b.U(), b.V(), gh);

    if ((blk.m + blk.n) * b.NcolUV() >= blk.m*blk.n) {
      // Better just to keep the whole block.
      b.U().Resize(0);
      b.V().Resize(0);
      if (ob && ob->HaveB())
        b.B() = ob->B();
      else if (gh)
        FillBFromGh(gh, CompressBlockInfo(blk.m, blk.n, opts.is_mrem,
                                          b.GetTol()),
                    blk, ma, b.B());
      else {
        Matrix<double> B(blk.m, blk.n);
        pma->Call(CompressBlockInfo(blk.m, blk.n, opts.is_mrem, b.GetTol()),
                  blk, NULL, NULL, B.GetPtr());
        b.B() = B;
      }
    } else if (opts.call_CompressQr) {
      double tol = (1.0 - opts.qr_alpha) * b.GetTol();
      CompressQr(b.U(), b.V(), !opts.is_mrem, tol, opts.allow_0rank);
    }
    if (ma_alloced) delete pma;
    if (gh) delete gh;
    return ret + 4;
  }}

  // Should never get here.
  assert(false);
  return -1;
}

//todo I think I should use a visitor pattern here to implement double
// dispatch on b and ob. For now, I'm going to commit the crime of using RTTI
// in the following two functions. I could implement single dispatch on b just
// by making a method wrapper to ApproxByLowRank for b, but as long as I'm
// commiting one crime, I'll just go ahead and commit two.
template<typename T>
inline int ApproxByLowRank (const LraOptions& opts, LraBlock* b,
                            MatrixAccessor& ma, const TypedLraBlock<T>* ob) {
  if (b->GetPrec() == 1)
    return ApproxByLowRank(opts, *(static_cast<TypedLraBlock<float>*>(b)),
                           ma, ob);
  else
    return ApproxByLowRank(opts, *(static_cast<TypedLraBlock<double>*>(b)),
                           ma, ob);
}

int ApproxByLowRank (const LraOptions& opts, LraBlock* b, MatrixAccessor& ma,
                     const LraBlock* ob)
  throw (OutOfMemoryException, UserReqException)
{
  if (ob) {
    if (ob->GetPrec() == 1)
      return ApproxByLowRank(opts, b, ma,
                             static_cast<const TypedLraBlock<float>*>(ob));
    else
      return ApproxByLowRank(opts, b, ma,
                             static_cast<const TypedLraBlock<double>*>(ob));
  } else {
    return ApproxByLowRank(opts, b, ma, (TypedLraBlock<float>*) NULL);
  }
}

// -----------------------------------------------------------------------------
// Compressor.

typedef unsigned long long int ull;

// Print to stdout "msg: 1 2 3 ... 100\n" to indicate progress.
class ProgressBar {
public:
  ProgressBar(const string& msg, ull total, ull output_lvl = 1);
  bool Incr(ull amount = 1);
  bool Finish();
  ull Progress() const { return _progress; }

private:
  ull _i, _total, _progress, _output_lvl;
};

ProgressBar::ProgressBar (const string& msg, ull total, ull output_lvl)
  : _i(0), _total(total), _progress(0), _output_lvl(output_lvl)
{
  if (_output_lvl > 0) {
    fprintf(stdout, "%s: ", msg.c_str());
    fflush(stdout);
  }
}

bool ProgressBar::Incr (ull amount) {
  if (_output_lvl > 0) {
    _i += amount;
    double p = 100.0 * (double) _i / _total;
    if ((ull) p > _progress) {
      _progress = (ull) p;
      // Don't print 100 here.
      if (_progress >= 100) return false;
      fprintf(stdout, "%lld ", _progress);
      fflush(stdout);
      return true;
    }
  }
  return false;
}

bool ProgressBar::Finish () {
  if (_output_lvl > 0) fprintf(stdout, "100\n");
  return true;
}

// Read an H-matrix sequentially one block at a time, return an LraBlock for
// each block as requested.
class HmatReader {
private:
  string _filename;
  FILE* _fid;
  double _tol;
  Blint _tol_denom;
  Blint _realp, _nb;
  bool _is_mrem, _is_compatible;
  vector<FileBlint> *_p, *_q;

public:
  HmatReader(const string& filename) throw (FileException);
  ~HmatReader();

  const string& GetFilename () const { return _filename; }
  double GetTol () const { return _tol; }

  // This checks everything.
  bool IsCompatible(Compressor::TolMethod tm, const vector<UInt>& p,
                    const vector<UInt>& q,
                    const vector<const MatBlock*>& sorted_blocks);
  LraBlock* ReadNextBlock() const;
};

HmatReader::HmatReader (const string& filename) throw (FileException)
  : _filename(filename), _is_compatible(false)
{
  _fid = fopen(filename.c_str(), "rb");
  if (!_fid) throw FileException("HmatReader");
  Blint m, n;
  _p = new vector<FileBlint>();
  _q = new vector<FileBlint>();
  Blint rerr_method;
  ReadHmatHeader(_fid, m, n, _realp, _nb, _tol, *_p, *_q, &rerr_method);
  _is_mrem = rerr_method == 1;
  if (!_is_mrem) // BREM
    _tol_denom = 1;
  else // MREM
    _tol_denom = m*n;
}

HmatReader::~HmatReader () {
  fclose(_fid);
  if (_p) delete _p;
  if (_q) delete _q;
}

// Comparison of MatBlock* based on r0, then c0.
static bool MatBlockLessThanRC (const MatBlock* b1, const MatBlock* b2) {
  if (b1->r0 < b2->r0) return true;
  if (b1->r0 > b2->r0) return false;
  return b1->c0 < b2->c0;
}

bool HmatReader::
IsCompatible (Compressor::TolMethod tm, const vector<UInt>& p,
              const vector<UInt>& q, const vector<const MatBlock*>& sbs) {
  if (_p) {
    // First check the header data.
    for (;;) { // breakable
      if (!(_is_compatible =
            (_is_mrem && (tm == Compressor::tm_mrem_fro ||
                          tm == Compressor::tm_mrem_abs)) ||
            (!_is_mrem && tm == Compressor::tm_brem_fro))) break;
      if (!(_is_compatible = p.size() == _p->size() &&
            q.size() == _q->size() && sbs.size() == (size_t) _nb)) break;
      for (size_t i = 0; i < p.size(); i++)
        if (!(_is_compatible = p[i] == (UInt) ((*_p)[i] + 1))) break;
      for (size_t i = 0; i < q.size(); i++)
        if (!(_is_compatible = q[i] == (UInt) ((*_q)[i] + 1))) break;
      break;
    }
    // We don't need these permutations anymore.
    delete _p; _p = NULL;
    delete _q; _q = NULL;
    if (!_is_compatible) return _is_compatible;

    // Checks so far, so now look at the blocks.
    for (UInt i = 0; i < sbs.size(); i++) {
      Blint r0, c0, m, n, rank;
      try {
        ReadHmatBlockInfo(_fid, _realp, r0, c0, m, n, rank);
      } catch (const Exception& e) {
        stringstream ss;
        ss << "IsCompatible: received [" << e.GetMsg() << "] on block " << i
           << " of " << sbs.size() << " blocks";
        throw FileException(ss.str());
      }
      MatBlock ob(r0, m, c0, n);
      vector<const MatBlock*>::const_iterator it =
        lower_bound(sbs.begin(), sbs.end(), &ob, MatBlockLessThanRC);
      const MatBlock* pb = *it;
      if (!(_is_compatible = (pb->r0 == (UInt) r0 && pb->c0 == (UInt) c0
                              && pb->m == (UInt) m && pb->n == (UInt) n)))
        break;
    }

    // Rewind to just after the header.
    fclose(_fid);
    _fid = fopen(_filename.c_str(), "rb");
    Blint m, n;
    ReadHmatHeader(_fid, m, n, _realp, _nb, _tol);
  }
  return _is_compatible;
}

template<typename T> static LraBlock*
ReadNextBlock (FILE* fid, bool is_mrem, double tol, Blint tol_denom) {
  Blint r0, c0, m, n, rank;
  T *B, *U, *Vt;
  ReadHmatBlock<T>(fid, r0, c0, m, n, rank, B, U, Vt);
  MatBlock b(r0, m, c0, n);
  if (is_mrem) tol *= sqrt((double) (m * n) / tol_denom);
  TypedLraBlock<T>* lb = new TypedLraBlock<T>(b, is_mrem, tol);
  if (B) {
    lb->B().Resize(m, n);
    memcpy(lb->B().GetPtr(), B, m*n*sizeof(T));
    delete[] B;
  } else {
    lb->U().Resize(m, rank);
    memcpy(lb->U().GetPtr(), U, m*rank*sizeof(T));
    delete[] U;
    lb->V().Resize(n, rank);
    Matrix<T> mVt(rank, n, Vt);
    Transpose(mVt, lb->V());
    delete[] Vt;
  }
  return lb;
}

LraBlock* HmatReader::ReadNextBlock () const
{
  if (_realp == 1)
    return hmmvp::ReadNextBlock<float>(_fid, _is_mrem, _tol, _tol_denom);
  else
    return hmmvp::ReadNextBlock<double>(_fid, _is_mrem, _tol, _tol_denom);
}

// Manager for the Parfor.
class CompressParforMgr : public mpi::ParforManager {
private:
  Compressor* _c;
  vector<mpi::ByteBufferWriter> _bw;
  ProgressBar _pb;
  vector<ull> _pb_incrs;

public:
  CompressParforMgr(Compressor* c);

  virtual int Isend(int job_idx, int dest);
  virtual int Irecv(int job_idx, int src) { return 0; }
  virtual void IsDone(int worker, int job_idx);
};

#define PB_ELEM 1

CompressParforMgr::CompressParforMgr (Compressor* c)
  : _c(c), _bw(mpi::GetNproc() - 1),
#if PB_ELEM
    _pb("Compress", (ull) _c->GetHmatNrows() * (ull) _c->GetHmatNcols(),
        c->GetOutputLevel()), _pb_incrs(mpi::GetNproc() - 1)
#else
    _pb("Compress", c->GetBlocks().size(), c->GetOutputLevel())
#endif
{}

int CompressParforMgr::Isend (int job_idx, int dest) {
  int idx = dest - 1;
  const MatBlock* b;
  LraBlock* lb = NULL;
  if (_c->HaveOldHmat()) {
    lb = _c->GetHmatReader()->ReadNextBlock();
    b = &lb->GetBlock();
  } else {
    b = _c->GetSortedBlocks()[job_idx];
  }
#if PB_ELEM
  _pb_incrs[idx] = (ull) b->m * (ull) b->n;
#endif
  _bw[idx].Reset();
  _bw[idx].WriteScalar<UInt>(b->r0);
  _bw[idx].WriteScalar<UInt>(b->m);
  _bw[idx].WriteScalar<UInt>(b->c0);
  _bw[idx].WriteScalar<UInt>(b->n);
  _bw[idx].WriteScalar<UInt>(lb ? lb->GetPrec() : 0);
  if (lb) {
    _bw[idx].WriteScalar((double) lb->GetTol());
    lb->Pack(_bw[idx]);
    delete lb;
  }
  _bw[idx].Isend(dest);
  return 0;
}

#if 0
#include <sys/resource.h>
#endif
static inline void ReportStuff (bool newline = true) {
#if 0
  static double initial = -1;
  rusage r;
  getrusage(RUSAGE_SELF, &r);
  if (initial < 0) { initial = r.ru_maxrss*1.0e-3; }
  errpr("%1.6f",
        //(r.ru_utime.tv_sec*1.0e6 + r.ru_utime.tv_usec)*1.0e-6,
        r.ru_maxrss*1.0e-3 - initial);
  if (newline) errpr("\n");
#endif
}

void CompressParforMgr::IsDone (int worker, int job_idx) {
  if (_pb.Incr(
#if PB_ELEM
               _pb_incrs[worker - 1]
#endif
               ))
    ReportStuff();
}

#undef PB_ELEM

// Worker for the Parfor.
class CompressParforWorker : public mpi::ParforWorker {
private:
  Compressor* _c;
  FILE* _fid;
  mpi::ByteBufferReaderMpi _br;

public:
  CompressParforWorker(Compressor* c, string fn) throw (FileException);
  virtual ~CompressParforWorker() { fclose(_fid); }

  virtual int Work(int job_idx, int root);
};

CompressParforWorker::CompressParforWorker (Compressor* c, string fn)
  throw (FileException)
  : _c(c)
{
  _fid = fopen(fn.c_str(), "wb");
  if (!_fid)
    throw FileException(string("CompressParforWorker: Could not open ") +
                        fn);
}

int CompressParforWorker::Work (int job_idx, int root) {
  MatBlock b;
  _br.Recv(root);
  b.r0 = _br.ViewScalarAndAdvance<UInt>();
  b.m  = _br.ViewScalarAndAdvance<UInt>();
  b.c0 = _br.ViewScalarAndAdvance<UInt>();
  b.n  = _br.ViewScalarAndAdvance<UInt>();
  UInt realp = _br.ViewScalarAndAdvance<UInt>();
  LraBlock* ob = NULL;
  if (realp) {
    double tol = _br.ViewScalarAndAdvance<double>();
    ob = NewLraBlock(b, realp, _c->IsMrem(), tol);
    ob->Unpack(_br);
  }
  _c->CompressBlockToFile(_fid, b, ob);
  if (ob) delete ob;
  return 0;
}

bool GreensFnMa::
Call (const CompressBlockInfo& cbi, const MatBlock& blk,
      const vector<UInt>* irs, const vector<UInt>* ics, double* B) const {
  vector<UInt> rs, cs;
  if (irs) {
    UInt n = irs->size();
    rs.resize(n);
    for (UInt i = 0; i < n; i++) rs[i] = _pp[(*irs)[i]];
  } else {
    rs.resize(blk.m);
    for (UInt i = 0; i < blk.m; i++) rs[i] = _pp[blk.r0 + i];
  }
  if (ics) {
    UInt n = ics->size();
    cs.resize(n);
    for (UInt i = 0; i < n; i++) cs[i] = _pq[(*ics)[i]];
  } else {
    cs.resize(blk.n);
    for (UInt i = 0; i < blk.n; i++) cs[i] = _pq[blk.c0 + i];
  }
#ifndef NDEBUG
  bool ret = _gf->Call(cbi, rs, cs, B);
  { size_t n = (irs ? irs->size() : blk.m) * (ics ? ics->size() : blk.n);
    for (size_t i = 0; i < n; i++) {
#if _MSC_VER > 1000
      assert(_finite(B[i]));
#else
      assert(!isnan(B[i]));
      assert(!isinf(B[i]));
#endif
    }}
  return ret;
#else
  return _gf->Call(cbi, rs, cs, B);
#endif
}

// We need to take care of the case of multiple Compressors running in
// parallel. We can handle this only in OpenMP. In MPI, only one Compressor may
// run at a time.

//todo-parallel In OpenMP, we are currently limited to multiple serial
// Compressors or one parallel one. With just a little work based around
// recognizing a tid other than 0 as root, we should be able to handle multiple
// parallel Compressors.

//todo-parallel In MPI, we can only handle one parallel Compressor. dc3dm's
// in-memory sub-Compressors give a lot of speedup, so until we provide that
// capability in MPI, 'dc3dm compress' must remain parallelized by OpenMP only.

/* This is a table of possible situations and the results. 'program' refers to
   how the overall program is running. 'Compressor' is how this particular
   Compressor is running. In the name 'AmRootOrIsSerial', 'Serial' refers to the
   Compressor. MPI and OpenMP cannot be built together at present, so it is not
   possible to have nt > 1 && np > 1.

     program  Compressor  (nt == 1 && np == 1) || id == 0  return
     serial   serial       t          t           t        t
     mpi      serial      not implemented
     openmp   serial
                root       t          t           t        t
            not root       t          t           f        t
     mpi         mpi
                root       t          f or t      t        t
            not root       t          f           f        f
     openmp   openmp
                root       f or t     t           t        t
            not root       f          t           f        f
*/
inline bool Compressor::AmRootOrIsSerial () const {
  return
    // Serial Compressor, serial Compressor in an OpenMP framework, or ...
    (_omp_nthreads == 1 && mpi::GetNproc() == 1) ||
    // ... anything, as long as I'm the root.
    mpi::AmTPid0();
}

void Compressor::GetHdData (const Hd* hd) {
  if (AmRootOrIsSerial()) {
    hd->Permutations(_ma.pp(), _ma.pq());
    if (mpi::GetNproc() > 1) {
      // Broadcast permutations.
      mpi::ByteBufferWriter bw;
      PackVec(bw, _ma.pp());
      PackVec(bw, _ma.pq());
      bw.Bcast();
    }
    // Only the root stores the blocks.
    _blocks.resize(hd->NbrBlocks());
    UInt i = 0;
    for (Hd::iterator it = hd->Begin(), end = hd->End(); it != end; ++it)
      _blocks[i++] = *it;
  } else {
    // Get broadcasted permutations.
    mpi::ByteBufferReaderMpi br;
    br.Bcast();
    UnpackVec(br, _ma.pp());
    UnpackVec(br, _ma.pq());
  }
}

static const vector<Blint>&
CopyAndDecr (const vector<UInt>& c, vector<Blint>& d) {
  d.resize(c.size());
  for (UInt i = 0; i < c.size(); i++) d[i] = c[i] - 1;
  return d;
}
    
void Compressor::WriteHmatHeader (FILE* fid) {
  FileBlint n;
  n = _ma.pp().size(); write(&n, 1, fid);
  n = _ma.pq().size(); write(&n, 1, fid);
  n = _prec; write(&n, 1, fid);
  // Switch from base-1 indexing (used by Hd) to base-0 indexing (used by
  // Hmat).
  vector<FileBlint> p;
  write(&CopyAndDecr(_ma.pp(), p)[0], _ma.pp().size(), fid);
  write(&CopyAndDecr(_ma.pq(), p)[0], _ma.pq().size(), fid);
  n = _blocks.size(); write(&n, 1, fid);
  // MREM (1) or BREM (0)?
  n = _tm == tm_brem_fro ? 0 : 1;
  write(&n, 1, fid);
  double d;
  d = _tol; write(&d, 1, fid);
  d = 14.0; write(&d, 1, fid);
}

Compressor::Compressor (const Hd* hd, GreensFn* gf)
  throw (Exception)
  : _ma(gf), _ohm(NULL), _have_old_hmat(false), _prec(2),
    _tm(Compressor::tm_mrem_abs), _Bfro(-1.0), _tol(1.0e-6),
    _user_set_prec(false), _avoid_redundant_gf_evals(false),
    _call_CompressQr(true), _allow_0rank(false), _output_lvl(1),
    _omp_nthreads(1), _ablr_stats(105, 0)
{
  if (!gf) throw Exception("gf is NULL.");
  GetHdData(hd);
}

Compressor::~Compressor () {
  if (_ohm) delete _ohm;
  if (AmRootOrIsSerial() && _output_lvl > 0)
    printf
      ("pid %2d:  "
       "old B %d  old QR %d  scratch B %d  scratch SVD %d  scratch ACA %d\n"
       "         B help SVD %d  B help ACA %d  UV help ACA %d\n"
       "         no help B %d  no help SVD %d  no help ACA %d\n",
       mpi::Pid(),
       _ablr_stats[1], _ablr_stats[2], _ablr_stats[0], _ablr_stats[3],
       _ablr_stats[4],
       _ablr_stats[13], _ablr_stats[24], _ablr_stats[34],
       _ablr_stats[100], _ablr_stats[103], _ablr_stats[104]);
}

Compressor* NewCompressor (const Hd* hd, GreensFn* gf) throw (Exception) {
  return new Compressor(hd, gf);
}

void DeleteCompressor (Compressor* c) { delete c; }

void Compressor::SetTolMethod (Compressor::TolMethod tm) { _tm = tm; }
Compressor::TolMethod Compressor::GetTolMethod () const { return _tm; }
bool Compressor::IsMrem () const {
  return _tm == tm_mrem_fro || _tm == tm_mrem_abs;
}

void Compressor::SetTol (double tol) throw (Exception) {
  if (tol <= 0.0) throw Exception("tol must be >= 0");
  _tol = tol;
}

void Compressor::SetPrec (UInt prec_code) throw (Exception) {
  if (!(prec_code == 1 || prec_code == 2))
    throw Exception("prec_code must be 1 or 2");
  _user_set_prec = true;
  _prec = prec_code;
}

void Compressor::Allow0RankBlocks (bool allow) { _allow_0rank = allow; }

void Compressor::SetOutputLevel (UInt lev) { _output_lvl = lev; }

void Compressor::AvoidRedundantGfCalls (bool flag) {
  _avoid_redundant_gf_evals = flag;
}

void Compressor::UseCompressQr (bool flag) { _call_CompressQr = flag; }

UInt Compressor::SetOmpNthreads (UInt n) {
  // We allow OpenMP only if there is just one process. It seems some
  // implementations don't play nicely together.
  if (n == 1) {
    _omp_nthreads = 1;
    return 1;
  }
  if (mpi::GetNproc() == 1) {
    assert(n > 0);
    omp_set_num_threads(n);
    _omp_nthreads = omp_get_max_threads();
  }
  return _omp_nthreads;
}

void Compressor::SetBfroEstimate (double Bfro)
  throw (Exception)
{
  if (Bfro <= 0.0) throw Exception("Bfro must be >= 0");
  _Bfro = Bfro;
}

double Compressor::GetBfroEstimate () const { return _Bfro; }

void Compressor::UseHmatFile (const string& hmat_filename)
  throw (Exception, FileException)
{
  bool can_read = true, is_compatible = true;
  if (AmRootOrIsSerial()) {
    try {
      _ohm = new HmatReader(hmat_filename);
    } catch (const FileException& e) {}
    if (!mpi_IsTrue(_ohm))
      can_read = false;
    else {
      // Sort the blocks for use in binary search.
      vector<const MatBlock*> sbs(_blocks.size());
      for (size_t i = 0; i < _blocks.size(); i++)
        sbs[i] = &_blocks[i];
      sort(sbs.begin(), sbs.end(), MatBlockLessThanRC);
      if (!mpi_IsTrue(_ohm->IsCompatible(_tm, _ma.pp(), _ma.pq(), sbs)))
        is_compatible = false;
      else {
        _have_old_hmat = true;
        mpi::Bcast(&_have_old_hmat, 1);
      }
    }
  } else {
    can_read = mpi_IsTrue();
    if (can_read) {
      is_compatible = mpi_IsTrue();
      if (is_compatible)
        mpi::Bcast(&_have_old_hmat, 1);
    }
  }
  if (!can_read)
    throw FileException("UseHmatFile: Could not read old H-matrix file.");
  if (!is_compatible)
    throw Exception("old H-matrix is not compatible with new one");
}

bool Compressor::HaveOldHmat () const { return _have_old_hmat; }

double Compressor::GetOldHmatBfro () const throw (Exception) {
  double Bfro = 0.0;
  if (!HaveOldHmat())
    throw Exception("GetOldHmatBfro: No old H-matrix.");
  if (AmRootOrIsSerial()) {
    Bfro = sqrt(HmatNormFrobenius2(_ohm->GetFilename()));
    if (_output_lvl > 0) printf("from old %e\n", Bfro);
    mpi::Bcast(&Bfro, 1);
  } else {
    mpi::Bcast(&Bfro, 1);
  }
  return Bfro;
}

static inline bool BlockIsOnDiag (const MatBlock& b) {
#if 0
  // This old definition works only for symmetric problems.
  return (b.r0 <= b.c0 + b.n - 1 && b.c0 <= b.r0 + b.m - 1);
#else
  // This method should work for all problems. It counts more than just
  // on-diagonal blocks, but as long as there aren't too many, this only
  // improves the accuracy of the estimate.
  return (b.m <= Hd::cluster_tree_min_points &&
          b.n <= Hd::cluster_tree_min_points);
#endif
}

static double FullBlockNormFro2 (MatrixAccessor& ma, const MatBlock& b) {
  Matrix<double> B(b.m, b.n);
  double* pB = B.GetPtr();
  //todo What's the correct tol for the cbi? For now, I'll just ask for exact
  // result, as this routine is used only for on-diagonal blocks.
  ma.Call(CompressBlockInfo(b.m, b.n, true, 0), b, NULL, NULL, pB);
  double nf2 = 0.0;
  for (UInt i = 0, mn = b.m*b.n; i < mn; i++) nf2 += pB[i]*pB[i];
  return nf2;
}

// Get a lower bound on ||B||_F based just on the diag blocks.
//todo Save these values and use them when compressing.
double Compressor::EstimateBfro ()
  throw (OutOfMemoryException, UserReqException)
{
  bool am_root = AmRootOrIsSerial();

  double nf2 = 0;
  if (mpi::GetNproc() == 1) { // serial or OpenMP
    ProgressBar pb("Estimate ||B||_F", _blocks.size(), GetOutputLevel());
    if (_omp_nthreads == 1) { // serial
      for (UInt i = 0, nb = _blocks.size(); i < nb; i++) {
        if (BlockIsOnDiag(_blocks[i]))
          nf2 += FullBlockNormFro2(_ma, _blocks[i]);
        pb.Incr();
      }
    } else { // OpenMP
      omp_set_num_threads(_omp_nthreads);
      int nthreads = omp_get_max_threads();
#pragma omp parallel
    { int tid = omp_get_thread_num();
      double nf2_me = 0.0;
      for (UInt i = tid, nb = _blocks.size(); i < nb; i += nthreads) {
        if (BlockIsOnDiag(_blocks[i]))
          nf2_me += FullBlockNormFro2(_ma, _blocks[i]);
        if (tid == 0) pb.Incr(nthreads);
      }
#pragma omp critical (nf2)
      nf2 += nf2_me; }
    }
    pb.Finish();
  } else { // parallel
    double nf2_me = 0;
    if (am_root) {
      // Assign blocks to procs.
      vector<UInt> idxs_me;
      UInt np = mpi::GetNproc();
      vector< mpi::ByteBufferWriter > _idxs(np - 1);
      for (UInt i = 0, k = 0, nb = _blocks.size(); i < nb; i++)
        if (BlockIsOnDiag(_blocks[i])) {
          if (k == (UInt) mpi::Root())
            idxs_me.push_back(i);
          else
            _idxs[k - 1].Write(&_blocks[i], 1);
          k = (k + 1) % np;
        }

      // Send out jobs.
      for (UInt i = 1; i < np; i++)
        _idxs[i - 1].Isend(i);

      // Do my part.
      ProgressBar pb("Estimate ||B||_F", idxs_me.size(), GetOutputLevel());
      for (size_t i = 0, n = idxs_me.size(); i < n; i++) {
        nf2_me += FullBlockNormFro2(_ma, _blocks[idxs_me[i]]);
        pb.Incr();
      }
      mpi::Barrier(); // Stay in scope until all msgs are delivered.
      pb.Finish();
    } else {
      // Receive my job.
      mpi::ByteBufferReaderMpi br;
      br.Recv(mpi::Root());
      UInt nb = br.Size() / sizeof(MatBlock);
      for (UInt i = 0; i < nb; i++) {
        MatBlock b;
        br.Read(&b, 1);
        nf2_me += FullBlockNormFro2(_ma, b);
      }
      mpi::Barrier(); // Match root's Barrier.
    }
    // Get the total.
    mpi::Allreduce(&nf2_me, &nf2, 1, MPI_SUM);
  }

  if (am_root && _output_lvl > 0) printf("est %e\n", sqrt(nf2));
  return sqrt(nf2);
}

static void GetTmpFilename (const string& base, string& fn) {
  stringstream ss;
  ss << base << "_" << mpi::Pid();
  fn = ss.str();
}

// Determine final abs tol and prec from _tol and _Bfro.
void Compressor::FinalizeTol () {
  if (_user_set_prec) return;
  // If the tolerance involves ||B||_F, include it.
  if (_tm == tm_mrem_fro) _tol *= _Bfro;
  if (IsMrem()) {
    // For each block, we need
    //   (tol Bfro)^2 blk.m blk.n / (M N) >= eps(_prec)^2 blk.m blk.n
    //    => (tol Bfro)^2 / (M N) >= eps(_prec)^2
    //    => _tol >= eps(_prec) sqrt(M N).
    // I use a factor of 10 to be conservative.
    _prec = (_tol < 10 * numeric_limits<float>::epsilon() *
             sqrt((double) _ma.pp().size() * _ma.pq().size())) ?
      GetPrecCode<double>() : GetPrecCode<float>();
  } else {
    // Block-wise relative error is easy. I use 10 again for padding.
    _prec = (_tol < 10 * numeric_limits<float>::epsilon()) ?
      GetPrecCode<double>() : GetPrecCode<float>();
  }
  if (_output_lvl > 0) cout << "prec " << _prec << endl;
}

// Send this many bytes at a time when concat'ing the files. I'm not sure
// whether one number is much better than any other.
const UInt Compressor::_send_sz = 1L << 24;

void Compressor::SendFile () {
  string fn;
  GetTmpFilename(_filename, fn);
  FILE* fid = fopen(fn.c_str(), "rb");

  mpi::ByteBufferWriter bw(sizeof(size_t) + _send_sz);
  for (;;) {
    bw.Reset();
    size_t nread =
      fread(bw.CopyTo<char>(_send_sz), sizeof(char), _send_sz, fid);
    bw.Resize<char>(nread);
    bw.Send(mpi::Root());
    if (nread < _send_sz * sizeof(char)) {
      // That was the last segment.
      if (nread > 0) {
        // Send a final empty message to signal I'm done.
        bw.Reset();
        bw.Send(mpi::Root());
      }
      break;
    }
  }

  fclose(fid);
}

void Compressor::
RecvAndWriteFile (FILE* fid, mpi::ByteBufferReaderMpi& br, UInt src) {
  for (;;) {
    br.Recv(src);
    if (br.Size() == 0) break;
    write(br.CurrentPosition(), br.Size(), fid);
  }
}

LraBlock* Compressor::
CompressBlock (const MatBlock& b, const LraBlock* ob, int& ret) {
  LraOptions opts;
  opts.is_mrem = IsMrem();
  double tol = _tol;
  if (opts.is_mrem) {
    // MREM formula for tol.
    tol *= sqrt((double) (b.m * b.n) / (_ma.pp().size() * _ma.pq().size()));
  }
  LraBlock* lb = NewLraBlock(b, _prec, opts.is_mrem, tol);
  if (_avoid_redundant_gf_evals) opts.min_aca_size = 2;
  opts.call_CompressQr = _call_CompressQr;
  opts.allow_0rank = _allow_0rank;
  _ma.StartBlock(CompressBlockInfo(b.m, b.n, opts.is_mrem, tol));
  if (ob)
    assert((ob->IsMrem() && _tm != tm_brem_fro) ||
           (!ob->IsMrem() && _tm == tm_brem_fro));
  ret = ApproxByLowRank(opts, lb, _ma, ob);
  _ma.FinishBlock(CompressBlockInfo(b.m, b.n, opts.is_mrem, tol));
  return lb;
}

void Compressor::
CompressBlockToFile (FILE* fid, const MatBlock& b, const LraBlock* ob) {
  int ret;
  LraBlock* lb = CompressBlock(b, ob, ret);
  if (_omp_nthreads == 1) {
    _ablr_stats[ret]++;
    lb->WriteToFile(fid);
  } else {
#pragma omp critical
    { _ablr_stats[ret]++;
      lb->WriteToFile(fid); }
  }
  delete lb;
}

void Compressor::
CompressBlockInMemory (HmatConstructor& hmc, size_t bi, const MatBlock& b,
                       const LraBlock* ob) {
  int ret;
  LraBlock* lb = CompressBlock(b, ob, ret);
  _ablr_stats[ret]++;
  lb->WriteInMemory(hmc, bi);
  delete lb;
}

// Serial or parallelized by OpenMP (as opposed to MPI).
void Compressor::CompressToFileSerial () {
  FILE* fid = fopen(_filename.c_str(), "ab");
  if (!fid) throw FileException(_filename + string(" can't be updated."));
  ProgressBar pb("Compress", _blocks.size(), GetOutputLevel());
  if (_have_old_hmat) {    
    for (UInt i = 0, nb = _blocks.size(); i < nb; i++) {
      LraBlock* ob = _ohm->ReadNextBlock();
      CompressBlockToFile(fid, ob->GetBlock(), ob);
      delete ob;
      pb.Incr();
    }
  } else {
    if (_omp_nthreads == 1) { // truly serial
      const size_t nb = _blocks.size();
      for (size_t i = 0; i < nb; i++) {
        CompressBlockToFile(fid, _blocks[i]);
        if (pb.Incr()) ReportStuff();
      }
    } else { // OpenMP
      omp_set_num_threads(_omp_nthreads);
      int nthreads = omp_get_max_threads();
      assert(nthreads > 0);
      SortBlocks();
      const size_t nb = _blocks.size();
#pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < nb; i++) {
        CompressBlockToFile(fid, *_sbs[i]);
        if (omp_get_thread_num() == 0 && pb.Incr(nthreads)) ReportStuff();
      }
    }
  }
  pb.Finish();
  fclose(fid);
}

void Compressor::ConcatFiles () {
  if (mpi::GetNproc() == 1) return;
  if (AmRootOrIsSerial()) {
    FILE* fid = fopen(_filename.c_str(), "ab");
    if (!fid) throw FileException(_filename + string(" can't be updated."));

    mpi::ByteBufferReaderMpi br;
    for (int pid = 1; pid < mpi::GetNproc(); pid++)
      RecvAndWriteFile(fid, br, pid);

    fclose(fid);
  } else {
    SendFile();
  }
}

// Comparison of MatBlock* based on m*n.
namespace {
  bool MatBlockGreaterThanMN (const MatBlock* b1, const MatBlock* b2)
  { return b1->m * b1->n > b2->m * b2->n; }
  bool MatBlockLessThanMN (const MatBlock* b1, const MatBlock* b2)
  { return b1->m * b1->n < b2->m * b2->n; }
}

void Compressor::SortBlocks () {
  _sbs.resize(_blocks.size());
  for (size_t i = 0; i < _blocks.size(); i++)
    _sbs[i] = &_blocks[i];
  // I go back and forth between this and MatBlockLessThanMN.
  sort(_sbs.begin(), _sbs.end(), MatBlockGreaterThanMN);
}

void Compressor::CompressToFileParallel () {
  CompressParforMgr* cpm = NULL;
  CompressParforWorker* cpw = NULL;
  bool all_ok;
  string msg = "";
  if (AmRootOrIsSerial()) {
    cpm = new CompressParforMgr(this);
    if (!HaveOldHmat()) SortBlocks();
    all_ok = mpi::AllOk(cpm);
  } else {
    string fn;
    GetTmpFilename(GetFilename(), fn);
    try {
      cpw = new CompressParforWorker(this, fn);
    } catch (const Exception& e) {
      stringstream ss;
      ss << "Pid " << mpi::Pid() << " got this when creating a worker: " <<
        e.GetMsg();
      msg = ss.str();
    }
    all_ok = mpi::AllOk(cpw);
  }
  if (all_ok) Parfor(cpm, cpw, _blocks.size());
  if (cpm) delete cpm;
  if (cpw) delete cpw;
  if (!all_ok) throw FileException(msg);
  ConcatFiles();    
}

void Compressor::CompressToFile (const string& hmat_fn)
  throw (OutOfMemoryException, UserReqException, FileException)
{
  _filename = hmat_fn;
  FinalizeTol();

  // Write the header.
  if (AmRootOrIsSerial()) {
    FILE* fid = fopen(_filename.c_str(), "wb");
    if (!fid) throw FileException(_filename + string(" can't be written."));
    WriteHmatHeader(fid);
    fclose(fid);
  }

  if (mpi::GetNproc() == 1) CompressToFileSerial();
  else CompressToFileParallel();
}

Hmat* Compressor::CompressInMemory (UInt ncol, UInt max_nthreads)
  throw (OutOfMemoryException, UserReqException)
{
  FinalizeTol();
  ProgressBar pb("Compress", _blocks.size(), GetOutputLevel());
  vector<Blint> p0(_ma.pp().size()), q0(_ma.pq().size());
  CopyAndDecr(_ma.pp(), p0);
  CopyAndDecr(_ma.pq(), q0);
  HmatConstructor hmc(_prec, _blocks.size(), p0, q0, ncol, max_nthreads);
  if (_have_old_hmat) {
    // It doesn't really make sense to compress an H-matrix in memory using an
    // old H-matrix file, but I'll keep this code for now.
    for (UInt i = 0, nb = _blocks.size(); i < nb; i++) {
      LraBlock* ob = _ohm->ReadNextBlock();
      CompressBlockInMemory(hmc, i, ob->GetBlock(), ob);
      delete ob;
      pb.Incr();
    }
  } else {
    if (_omp_nthreads == 1) { // truly serial
      const size_t nb = _blocks.size();
      for (size_t i = 0; i < nb; i++) {
        CompressBlockInMemory(hmc, i, _blocks[i]);
        if (pb.Incr()) ReportStuff();
      }
    } else { // OpenMP
      omp_set_num_threads(_omp_nthreads);
      int nthreads = omp_get_max_threads();
      assert(nthreads > 0);
      SortBlocks();
      const size_t nb = _blocks.size();
#pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < nb; i++) {
        CompressBlockInMemory(hmc, i, *_sbs[i]);
        if (AmRootOrIsSerial() && pb.Incr(nthreads)) ReportStuff();
      }
    }
    pb.Finish();
  }
  return hmc.GetHmat();
}
}
