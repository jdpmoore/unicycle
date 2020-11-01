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

#include <assert.h>
#include <math.h>
#include <list>
#include <limits>
#include <algorithm>
#include "util/include/LinAlg.hpp"

namespace hmmvp {
using namespace util;

// -----------------------------------------------------------------------------
// Recompress using skinny QR (also sometimes known as low-rank SVD, I think).

template<typename T>
int SvdChooseRank (const Matrix<T>& s, bool use_rel_err, double err,
                   bool allow_0rank) {
  int n = s.Size();
  double fac = 1;
  const T* rs = s.GetPtr();
  if (use_rel_err) {
    fac = 0;
    for (int i = 0; i < n; i++) fac += rs[i]*rs[i];
  }
  int rank;
  double cs = 0, err2 = err*err;
  for (rank = n; rank > 0; rank--) {
    T se = rs[rank - 1];
    cs += se*se;
    if (cs > err2*fac) break;
  }
  // In practice, I prefer not to have rank-0 blocks.
  if (!allow_0rank && rank == 0 && n > 0) rank = 1;
  return rank;
}

template<typename T>
void CompressQr (Matrix<T>& U, Matrix<T>& V, const bool use_rel_err,
                 const double err, bool allow_0rank) {
  assert (U.Size(2) == V.Size(2));
  int n = U.Size(2);
  if (n == 0) return;
  int rank;
  Matrix<T> Uq, Vq, U1, s, V1t;
  { Matrix<T> C(n, n);
    { Matrix<T> Ur(n, n), Vr(n, n);
      SkinnyQr(U, Uq, Ur);
      SkinnyQr(V, Vq, Vr);
      Mult(Ur, false, Vr, true, C); }
    Matrix<T> V1to;
    Svd(C, U1, s, V1to);
    rank = SvdChooseRank(s, use_rel_err, err, allow_0rank);
    if (rank == n) return;
    U1.Reshape(n, rank);
    V1t.Resize(rank, n);
    T* pV1t = V1t.GetPtr();
    T* pV1to = V1to.GetPtr();
    for (int j = 0; j < n; j++) {
      memcpy(pV1t, pV1to, rank*sizeof(T));
      pV1t += rank;
      pV1to += n;
    }
  }
  if (rank == 0) return;

  // U1(:,1:rank)*diag(s(1:rank))
  T* dU1 = U1.GetPtr();
  T* ds = s.GetPtr();
  for (int j = 0; j < rank; j++) {
    T sj = ds[j];
    for (int i = 0; i < n; i++) dU1[i] *= sj;
    dU1 += n;
  }

  Mult(Uq, false, U1, false, U);
  Mult(Vq, false, V1t, true, V);
}

// -----------------------------------------------------------------------------
// Implementation of Bebendorf's ACA.

// x = U(:,1:n-1)' U(:,n). n >= 2 is required.
template<typename T>
static void CalcUtu (const list< Matrix<T> >& U, Matrix<T>& x) {
  const Matrix<T>& u = U.back();
  int m = u.Size(1);
  x.Resize(U.size() - 1);
  int i = 1;
  typename list< Matrix<T> >::const_iterator
    it = U.begin(),
    end = U.end();
  --end;
  for (; it != end; ++it, i++)
    x(i) = dot(m, it->GetPtr(), 1, u.GetPtr(), 1);
}

template<typename T>
static inline double NormFro2uv (const Matrix<T>& u, const Matrix<T>& v) {
  int Um = u.Size(1), Vm = v.Size(1);
  return dot(Um, u.GetPtr(), 1, u.GetPtr(), 1) *
         dot(Vm, v.GetPtr(), 1, v.GetPtr(), 1);
}

// Calculate the contribution of U(:,end), V(:,end) to norm(U V', 'fro')^2.
template<typename T> static void
UpdateNormFro2UV (const list< Matrix<T> >& U, const list< Matrix<T> >& V,
                  double& nf2, double& u2v2) {
  int n = U.size();
  if (n == 0) return;
  const Matrix<T>& u = U.back();
  const Matrix<T>& v = V.back();
  u2v2 = NormFro2uv(u, v);
  nf2 += u2v2;
  if (n == 1) return;
  Matrix<T> Utu, Vtv;
  CalcUtu(U, Utu);
  CalcUtu(V, Vtv);
  nf2 += 2*dot(n - 1, Utu.GetPtr(), 1, Vtv.GetPtr(), 1);
}

template<typename T>
static inline Matrix<T>& ListMatPushBackZeroVec (list< Matrix<T> >& lm, int m) {
  lm.push_back(Matrix<T>());
  lm.back().Reshape(m);
  return lm.back();
}

template<typename T>
static inline void ListMatPushBackVec 
(list< Matrix<T> >& lm, double* v, const int m) {
  lm.push_back(Matrix<T>());
  Matrix<T>& vm = lm.back();
  vm.Resize(m);
  T* pvm = vm.GetPtr();
  for (int i = 0; i < m; i++) pvm[i] = (T) v[i];
}

// The list is emptied as the matrix is formed. Optionally transpose the
// matrix. Each matrix in lv must have the same size.
template<typename T1, typename T2>
static void ListMatToMat (list< Matrix<T1> >& lv, bool want_transp,
                          Matrix<T2>& A) {
  int n = lv.size();
  if (n == 0) {
    A.Resize(0);
    return;
  }
  int m = lv.begin()->Size();
  if (want_transp) {
    A.Resize(n, m);
    T2* pA = A.GetPtr();
    int i = 0;
    for (typename list< Matrix<T1> >::iterator it = lv.begin(), end = lv.end();
         it != end; ++it) {
      T1* pc = it->GetPtr();
      for (int j = 0; j < m; j++) pA[j*n + i] = (T2) pc[j];
      i++;
      it->Resize(0);
    }
  } else {
    A.Resize(m, n);
    T2* pA = A.GetPtr();
    int i = 0;
    for (typename list< Matrix<T1> >::iterator it = lv.begin(), end = lv.end();
         it != end; ++it) {
      T1* pc = it->GetPtr();
      for (int j = 0; j < m; j++) pA[j] = (T2) pc[j];
      i++;
      it->Resize(0);
      pA += m;
    }
  }
}

template<typename T>
static double ScaleFactor () {
  double scale_fac = numeric_limits<T>::epsilon();
  if (sizeof(T) == sizeof(float)) return 1.0e1 * scale_fac;
  else return 1.0e2 * scale_fac;
}

template<typename T>
static int GetMaxElement (const T* v, const int n) {
  T ma = 0;
  int im = 0;
  for (int i = 0; i < n; i++) {
    T av = fabs(v[i]);
    if (av > ma) {
      im = i;
      ma = av;
    }
  }
  return im;
}

template<typename T>
static int GetInitialIrow (const Matrix<T>& U, UInt m) {
  int n = U.Size(2);
  int irow;
  if (n == 0)
    irow = m/2;
  else {
    int m = U.Size(1);
    // Find the max element in the last column.
    irow = GetMaxElement(U.GetPtr() + m*(n - 1), m);
  }
  return irow;
}

// Fills lU, lV. Returns norm(U V', 'fro')^2.
template<typename T1, typename T2>
static double GetInitialUV (const Matrix<T2>& U, const Matrix<T2>& V,
                            list< Matrix<T1> >& lU, list< Matrix<T1> >& lV) {
  int Um = U.Size(1), Vm = V.Size(1), n = U.Size(2);
  const T2* pU = U.GetPtr();
  const T2* pV = V.GetPtr();
  double nf2 = 0;
  for (int i = 1; i <= n; i++) {
    T1* u = ListMatPushBackZeroVec(lU, Um).GetPtr();
    T1* v = ListMatPushBackZeroVec(lV, Vm).GetPtr();
    for (int j = 0; j < Um; j++) u[j] = (T1) pU[j];
    for (int j = 0; j < Vm; j++) v[j] = (T1) pV[j];
    pU += Um;
    pV += Vm;
    double u2v2;
    UpdateNormFro2UV(lU, lV, nf2, u2v2);
  }
  return nf2;
}

// a -= A*B(ib+1,:).'
template<typename T>
static void CalcamAb (double* a, const list< Matrix<T> >& A,
                      const list< Matrix<T> >& B, int ib, double* wrk) {
  bool first = true;
  UInt m = 0;
  for (typename list< Matrix<T> >::const_iterator
         Ait = A.begin(), Aend = A.end(), Bit = B.begin();
       Ait != Aend; ++Ait, ++Bit) {
    if (first) {
      first = false;
      m = Ait->Size();
      memset(wrk, 0, m*sizeof(double));
    }
    double b = Bit->GetPtr()[ib];
    const T* Ac = Ait->GetPtr();
    for (UInt i = 0; i < m; i++) wrk[i] += ((double) Ac[i])*b;
  }
  if (!first) {
    for (UInt i = 0; i < m; i++) a[i] -= wrk[i];
  }
}

static inline int
MaxAvailable (double* u, const int m, const vector<char>& unused_rows)
throw (Exception) {
  assert (m > 0);
  double ma = 0;
  int irow = -1;
  for (int i = 0; i < m; i++)
    if (unused_rows[i]) {
      double au = fabs(u[i]);
      if (au >= ma) {
        ma = au;
        irow = i;
      }
    }
  if (irow < 0) throw Exception();
  return irow;
}

class AcaGfHolder {
private:
  const MatBlock* _blk;
  vector<double*> _rows, _cols;
  vector<double> _dwrk;
  vector<UInt> _uiwrk;

  void FreeVecDoublePtr (vector<double*>& v) {
    for (size_t i = 0; i < v.size(); i++) if (v[i]) delete[] v[i];
  }
  //todo-opt Could avoid repeatedly checking and indexing.
  bool ValueAt (UInt r, UInt c, double* d) {
    if (_rows[r]) { *d = _rows[r][c]; return true; }
    if (_cols[c]) { *d = _cols[c][r]; return true; }
    return false;
  }
  void InsertRow (UInt r, const double* data) {
    assert(!_rows[r]);
    _rows[r] = new double[_blk->n];
    memcpy(_rows[r], data, _blk->n * sizeof(*data));
  }
  void InsertCol (UInt c, const double* data) {
    // Originally I had this assert. However, it can't actually always be true
    // because numerical errors in ACA can make a column be requested multiple
    // times; for, unlike rows, logic does not control access to the columns,
    // but only finite precision arithmetic.
    //assert(!_cols[c]);
    // So instead I check.
    if (!_cols[c]) {
      _cols[c] = new double[_blk->m];
      memcpy(_cols[c], data, _blk->m * sizeof(*data));
    }
  }

public:
  AcaGfHolder (const MatBlock* blk)
    : _blk(blk), _rows(blk->m, NULL), _cols(blk->n, NULL) {
    UInt K = max(blk->m, blk->n);
    _dwrk.resize(K);
    _uiwrk.reserve(K);
  }
  ~AcaGfHolder () {
    FreeVecDoublePtr(_rows);
    FreeVecDoublePtr(_cols);
  }

  bool Call (const CompressBlockInfo& cbi, MatrixAccessor& ma,
             const vector<UInt>* row, const vector<UInt>* col, double *data,
             bool insert = true) {
    _uiwrk.clear();
    if (row) {
      for (size_t i = 0; i < _blk->n; i++)
        if (!ValueAt((*row)[0] - _blk->r0, i, data + i))
          _uiwrk.push_back(_blk->c0 + i);
      ma.Call(cbi, *_blk, row, &_uiwrk, &_dwrk[0]);
      for (size_t i = 0; i < _uiwrk.size(); i++)
        data[_uiwrk[i] - _blk->c0] = _dwrk[i];
      if (insert) InsertRow((*row)[0] - _blk->r0, data);
    } else {
      for (size_t i = 0; i < _blk->m; i++)
        if (!ValueAt(i, (*col)[0] - _blk->c0, data + i))
          _uiwrk.push_back(_blk->r0 + i);
      ma.Call(cbi, *_blk, &_uiwrk, col, &_dwrk[0]);
      for (size_t i = 0; i < _uiwrk.size(); i++)
        data[_uiwrk[i] - _blk->r0] = _dwrk[i];
      if (insert) InsertCol((*col)[0] - _blk->c0, data);
    }
#if 0
    ma.Call(cbi, *_blk, row, col, &_dwrk[0]);
    double err = 0;
    for (size_t i = 0, n = row ? _blk->n : _blk->m; i < n; i++) {
      double diff = data[i] - _dwrk[i];
      err += diff*diff;
    }
    if (err != 0) printf("m %ld n %ld err %1.2e\n", _blk->m, _blk->n, err);
#endif
    return true;
  }

  //dbg
  bool HaveCol (size_t c) { return _cols[c]; }
};

// Essentially a duplicate of ACAr in Bebendorf's AHMED/ACA.h.
template<typename T>
void Aca (const MatBlock& blk, MatrixAccessor& ma, double scale,
          const bool use_rel_err, const double err,
          Matrix<T>& U, Matrix<T>& V, AcaGfHolder* gh)
  throw (OutOfMemoryException, UserReqException)
{
  CompressBlockInfo cbi(blk.m, blk.n, !use_rel_err, err);
  double scale_fac = ScaleFactor<T>();
  int k1 = U.Size(2);
  int irow = GetInitialIrow(U, blk.m);

  // Use double for both T = float and double.
  list< Matrix<double> > lU, lV;
  double Bfro2 = GetInitialUV(U, V, lU, lV);
  U.Resize(0);
  V.Resize(0);
  // If this is a warm start, come up with a better scale.
  if (Bfro2 > 0) scale = sqrt(Bfro2/(blk.m*blk.n));

  // ma call results go here. Set up some aliases.
  UInt mamn = std::max(blk.m, blk.n);
  vector<double> rwrk(2*mamn);
  double *ru = &rwrk[0], *rv = ru, *wrk = ru + mamn;

  vector<char> unused_rows(blk.m, 1);
  vector<UInt> si(1);
  for (int k = k1, k_max = std::min(blk.m, blk.n); k < k_max; k++) {
    // Get a new row of the remainder matrix R.
    bool fnd_irow = false;
    int icol;
    for (;;) {
      si[0] = blk.r0 + irow;
      bool success = gh ?
        gh->Call(cbi, ma, &si, NULL, rv) :
        ma.Call(cbi, blk, &si, NULL, rv);
      if (!success) {
        if (gh) delete gh;
        throw UserReqException("in ACA");
      }
      // v -= U(irow,:)*V.'
      CalcamAb(rv, lV, lU, irow, wrk);
      // New column.
      icol = GetMaxElement(rv, blk.n);
      unused_rows[irow] = 0;
      // Acceptable?
      fnd_irow = fabs(rv[icol]) > scale_fac * scale;
      if (fnd_irow) break;
      // No, so try a new irow.
      bool fnd = false;
      for (UInt i = 0; i < blk.m; i++) {
        irow = (irow + 1) % blk.m;
        if (unused_rows[irow]) { fnd = true; break; }
      }
      if (!fnd) break;
      // Increment k because we have one fewer slot in unused_rows.
      k++;
    }
    // Check for exact 0.
    if (fnd_irow) fnd_irow = rv[icol] != 0;
    // Exactly (to measured precision) low-rank matrix.
    if (!fnd_irow) break;
    double den = rv[icol];
    for (UInt i = 0; i < blk.n; i++) rv[i] /= den;
    ListMatPushBackVec(lV, rv, blk.n);
    // Get the column of R.
    si[0] = blk.c0 + icol;
    bool success = gh ?
      gh->Call(cbi, ma, NULL, &si, ru) :
      ma.Call(cbi, blk, NULL, &si, ru);
    if (!success) throw UserReqException("in ACA");
    // u -= U*V(icol,:).'
    CalcamAb(ru, lU, lV, icol, wrk);
    ListMatPushBackVec(lU, ru, blk.m);
    if (k >= k_max - 1) break;
    // Next irow.
    try {
      irow = MaxAvailable(ru, blk.m, unused_rows);
    } catch (const Exception& e) {
      double mu = 0;
      for (UInt i = 0; i < blk.m; i++) mu = std::max(mu, fabs(ru[i]));
      errpr("k %d m %ld n %ld scale %e scale_fac %e irow %d icol %d "
            "max |ru| %e err %e\n",
            k, blk.m, blk.n, scale, scale_fac, irow, icol, mu, err);
      if (blk.m * blk.n <= 10000) {
        Matrix<double> G1(blk.m, blk.n);
        vector<UInt> rs(blk.m), cs(blk.n);
        for (size_t i = 0; i < blk.m; i++) rs[i] = blk.r0 + i;
        for (size_t i = 0; i < blk.n; i++) cs[i] = blk.c0 + i;
        ma.Call(cbi, blk, &rs, &cs, G1.GetPtr());
        cout << "G1 " << G1 << endl;
      }
      throw UserReqException("Aca: irow exception");
    }
    // Assess error.
    double Rfro2;
    if (use_rel_err) {
      UpdateNormFro2UV(lU, lV, Bfro2, Rfro2);
      if (Rfro2 <= err*err*Bfro2) break;
    } else {
      Rfro2 = NormFro2uv(lU.back(), lV.back());
      if (Rfro2 <= err*err) break;
    }
    scale = sqrt(Rfro2 / (blk.m*blk.n));
  }

  ListMatToMat(lU, false, U);
  ListMatToMat(lV, false, V);
}

// ---------------------------------------------------------------------------

template<typename T>
TypedLraBlock<T>::TypedLraBlock (const MatBlock& b, bool is_mrem, double tol)
  : LraBlock(b, is_mrem, tol) {}

template<> inline UInt TypedLraBlock<float >::GetPrec () const { return 1; }
template<> inline UInt TypedLraBlock<double>::GetPrec () const { return 2; }

}
