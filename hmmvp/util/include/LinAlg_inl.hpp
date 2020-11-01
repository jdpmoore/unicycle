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

namespace util {
extern "C" void sgemm_(
  char*, char*, blas_int*, blas_int*, blas_int*, float*, const float*,
  blas_int*, const float*, blas_int*, float*, float*, blas_int*);
extern "C" void sgemv_(
  char*, blas_int*, blas_int*, float*, const float*, blas_int*, const float*,
  blas_int*, float*, float*, blas_int*);
extern "C" float sdot_(
  blas_int*, const float*, blas_int*, const float*, blas_int*);
extern "C" void saxpy_(
  blas_int*, float*, const float*, blas_int*, float*, blas_int*);
extern "C" void sgesvd_(
  char* jobu, char* jobvt, blas_int* m, blas_int* n, float* A, blas_int* lda,
  float* s, float* U, blas_int* ldu, float* Vt, blas_int* ldvt, float* wrk,
  blas_int* lwrk, blas_int* info);
extern "C" void sgeqrf_(
  blas_int* m, blas_int* n, float* A, blas_int* lda, float* tau, float* work,
  blas_int* lwork, blas_int* info);
extern "C" void sorgqr_(
  blas_int* m, blas_int* n, blas_int* k, float* A, blas_int* lda, float* tau,
  float* work, blas_int* lwork, blas_int* info);
extern "C" void sgelqf_(
  blas_int* m, blas_int* n, float* A, blas_int* lda, float* tau, float* work,
  blas_int* lwork, blas_int* info);
extern "C" void sorglq_(
  blas_int* m, blas_int* n, blas_int* k, float* A, blas_int* lda, float* tau,
  float* work, blas_int* lwork, blas_int* info);
extern "C" void strtrs_(
  char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, float* A,
  blas_int* lda, float* B, blas_int* ldb, blas_int* info);

extern "C" void dgemm_(
  char*, char*, blas_int*, blas_int*, blas_int*, double*, const double*,
  blas_int*, const double*, blas_int*, double*, double*, blas_int*);
extern "C" void dgemv_(
  char*, blas_int*, blas_int*, double*, const double*, blas_int*,
  const double*, blas_int*, double*, double*, blas_int*);
extern "C" double ddot_(
  blas_int*, const double*, blas_int*, const double*, blas_int*);
extern "C" void daxpy_(
  blas_int*, double*, const double*, blas_int*, double*, blas_int*);
extern "C" void dgesvd_(
  char* jobu, char* jobvt, blas_int* m, blas_int* n, double* A, blas_int* lda,
  double* s, double* U, blas_int* ldu, double* Vt, blas_int* ldvt, double* wrk,
  blas_int* lwrk, blas_int* info);
extern "C" void dgeqrf_(
  blas_int* m, blas_int* n, double* A, blas_int* lda, double* tau,
  double* work, blas_int* lwork, blas_int* info);
extern "C" void dorgqr_(
  blas_int* m, blas_int* n, blas_int* k, double* A, blas_int* lda, double* tau,
  double* work, blas_int* lwork, blas_int* info);
extern "C" void dgelqf_(
  blas_int* m, blas_int* n, double* A, blas_int* lda, double* tau,
  double* work, blas_int* lwork, blas_int* info);
extern "C" void dorglq_(
  blas_int* m, blas_int* n, blas_int* k, double* A, blas_int* lda, double* tau,
  double* work, blas_int* lwork, blas_int* info);
extern "C" void dtrtrs_(
  char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, double* A,
  blas_int* lda, double* B, blas_int* ldb, blas_int* info);
  
template<> inline void gemm<double> (
  char transa, char transb, blas_int m, blas_int n, blas_int k, double alpha,
  const double* A, blas_int lda, const double* B, blas_int ldb, double beta,
  double* C, blas_int ldc)
{
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
         C, &ldc);
}

template<> inline void gemm<float> (
  char transa, char transb, blas_int m, blas_int n, blas_int k,
  float alpha, const float* A, blas_int lda,
  const float* B, blas_int ldb, float beta,
  float* C, blas_int ldc)
{
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
         C, &ldc);
}

template<> inline void gemv<double> (
  char trans, blas_int m, blas_int n, double alpha,
  const double* a, blas_int lda, const double* x, blas_int incx, double beta,
  double* y, blas_int incy)
{
  dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<> inline void gemv<float> (
  char trans, blas_int m, blas_int n, float alpha,
  const float* a, blas_int lda, const float* x, blas_int incx, float beta,
  float* y, blas_int incy)
{
  sgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<> inline double dot<double> (
  blas_int n, const double* x, blas_int incx, const double* y, blas_int incy)
{
  return ddot_(&n, x, &incx, y, &incy);
}

template<> inline float dot<float> (
  blas_int n, const float* x, blas_int incx, const float* y, blas_int incy)
{
  return sdot_(&n, x, &incx, y, &incy);
}

template<> inline void axpy<double> (
  blas_int n, double a, const double* x, blas_int incx, double* y,
  blas_int incy)
{
  daxpy_(&n, &a, x, &incx, y, &incy);
}

template<> inline void axpy<float> (
  blas_int n, float a, const float* x, blas_int incx, float* y, blas_int incy)
{
  saxpy_(&n, &a, x, &incx, y, &incy);
}

template<> inline void gesvd<float> (
  char jobu, char jobvt, blas_int m, blas_int n, float* A, blas_int lda,
  float* s, float* U, blas_int ldu, float* Vt, blas_int ldvt, float* wrk,
  blas_int lwrk, blas_int& info)
{
  info = 0;
  sgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, Vt, &ldvt, wrk, &lwrk,
          &info);
}

template<> inline void gesvd<double> (
  char jobu, char jobvt, blas_int m, blas_int n, double* A, blas_int lda,
  double* s, double* U, blas_int ldu, double* Vt, blas_int ldvt, double* wrk,
  blas_int lwrk, blas_int& info)
{
  info = 0;
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, Vt, &ldvt, wrk, &lwrk,
          &info);
}

template<> inline void geqrf<float> (
  blas_int m, blas_int n, float* A, blas_int lda, float* tau, float* work,
  blas_int lwork, blas_int& info)
{
  info = 0;
  sgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
}

template<> inline void geqrf<double> (
  blas_int m, blas_int n, double* A, blas_int lda, double* tau, double* work,
  blas_int lwork, blas_int& info)
{
  info = 0;
  dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
}

template<> inline void gelqf<float> (
  blas_int m, blas_int n, float* A, blas_int lda, float* tau, float* work,
  blas_int lwork, blas_int& info)
{
  info = 0;
  sgelqf_(&m, &n, A, &lda, tau, work, &lwork, &info);
}

template<> inline void gelqf<double> (
  blas_int m, blas_int n, double* A, blas_int lda, double* tau, double* work,
  blas_int lwork, blas_int& info)
{
  info = 0;
  dgelqf_(&m, &n, A, &lda, tau, work, &lwork, &info);
}

template<> inline void orgqr<float> (
  blas_int m, blas_int n, blas_int k, float* A, blas_int lda, float* tau,
  float* work, blas_int lwork, blas_int& info)
{
  info = 0;
  sorgqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
}

template<> inline void orgqr<double> (
  blas_int m, blas_int n, blas_int k, double* A, blas_int lda, double* tau,
  double* work, blas_int lwork, blas_int& info)
{
  info = 0;
  dorgqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
}

template<> inline void orglq<float> (
  blas_int m, blas_int n, blas_int k, float* A, blas_int lda, float* tau,
  float* work, blas_int lwork, blas_int& info)
{
  info = 0;
  sorglq_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
}

template<> inline void orglq<double> (
  blas_int m, blas_int n, blas_int k, double* A, blas_int lda, double* tau,
  double* work, blas_int lwork, blas_int& info)
{
  info = 0;
  dorglq_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
}

template<> inline void trtrs<float> (
  char uplo, char trans, char diag, blas_int n, blas_int nrhs, float* A,
  blas_int lda, float* B, blas_int ldb, blas_int& info)
{
  info = 0;
  strtrs_(&uplo, &trans, &diag, &n, &nrhs, A, &lda, B, &ldb, &info);
}

template<> inline void trtrs<double> (
  char uplo, char trans, char diag, blas_int n, blas_int nrhs, double* A,
  blas_int lda, double* B, blas_int ldb, blas_int& info)
{
  info = 0;
  dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, A, &lda, B, &ldb, &info);
}

template<typename T>
void Mult (const Matrix<T>& A, bool do_trans_A, const Matrix<T>& B,
           bool do_trans_B, Matrix<T>& C) {
  blas_int m, n, k;
  char transa = 'n', transb = 'n';
  if (do_trans_A) {
    transa = 't';
    m = A.Size(2);
    k = A.Size(1);
  } else {
    m = A.Size(1);
    k = A.Size(2);
  }
  if (do_trans_B) {
    transb = 't';
    n = B.Size(1);
  } else {
    n = B.Size(2);
  }
  blas_int lda = A.Size(1), ldb = B.Size(1);
  C.Resize(m, n);
  gemm(transa, transb, m, n, k, (T) 1.0, A.GetPtr(), lda, B.GetPtr(), ldb,
       (T) 0.0, C.GetPtr(), m);
}

template<typename T1, typename T2>
void Svd (const Matrix<T1>& A, Matrix<T2>& U, Matrix<T2>& s, Matrix<T2>& Vt)
  throw (LapackException)
{
  if (A.Size() == 0) throw LapackException("Svd: A.Size() == 0");
  blas_int m, lda, n, ldu, ldvt, info, lwork;
  m = lda = ldu = A.Size(1);
  n = A.Size(2);
  ldvt = std::min(m, n);
  U.Resize(m, ldvt);
  s.Resize(ldvt);
  Vt.Resize(ldvt, n);
  char jobu = 'S', jobvt = 'S';
  lwork = -1; // workspace query
  T2 dlwork;
  gesvd(jobu, jobvt, m, n, (T2*) NULL, lda, s.GetPtr(), U.GetPtr(), ldu,
        Vt.GetPtr(), ldvt, &dlwork, lwork, info);
  // A copy, U, s, Vt, work
  lwork = (blas_int) dlwork;
  WorkArray<T2> w(lwork + m*n);
  Matrix<T2> work;
  work.SetPtr(lwork, w.AllocWork(lwork));
  Matrix<T2> A1;
  A1.SetPtr(m, n, w.AllocWork(m*n));
  A1 = A;
  gesvd(jobu, jobvt, m, n, A1.GetPtr(), lda, s.GetPtr(), U.GetPtr(), ldu,
        Vt.GetPtr(), ldvt, work.GetPtr(), lwork, info);
}

template<typename T>
void SkinnyQr (Matrix<T>& A, Matrix<T>& R) throw (LapackException) {
  if (A.Size() == 0) throw LapackException("SkinyyQr: A.Size() == 0");
  blas_int m, lda, n, k, info, lwork;
  m = lda = A.Size(1);
  n = A.Size(2);
  if (m < n) throw LapackException("SkinnyQr: m < n");
  k = n;
  R.Resize(k, k);
  lwork = -1;
  T dlwork;
  geqrf(m, n, (T*) NULL, lda, (T*) NULL, &dlwork, lwork, info);
  blas_int lwork1 = (blas_int) dlwork;
  orgqr(m, n, k, (T*) NULL, lda, (T*) NULL, &dlwork, lwork, info);
  blas_int lwork2 = (blas_int) dlwork;
  lwork = std::max(lwork1, lwork2);
  WorkArray<T> w(lwork + k);
  Matrix<T> work, tau;
  work.SetPtr(lwork, w.AllocWork(lwork));
  tau.SetPtr(k, w.AllocWork(k));
  geqrf(m, n, A.GetPtr(), lda, tau.GetPtr(), work.GetPtr(), lwork, info);
  T* rQ = A.GetPtr();
  T* rR = R.GetPtr();
  memset(rR, 0, k*k*sizeof(T));
  for (int c = 0; c < k; c++) {
    memcpy(rR, rQ, (c + 1)*sizeof(T));
    rR += n;
    rQ += m;
  }
  orgqr(m, n, k, A.GetPtr(), lda, tau.GetPtr(), work.GetPtr(), lwork, info);
}

template<typename T>
void SkinnyQr (const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R)
  throw (LapackException)
{
  Q = A;
  SkinnyQr(Q, R);
}

template<typename T>
void ShortLq (Matrix<T>& A, Matrix<T>& L) throw (LapackException) {
  if (A.Size() == 0) throw LapackException("SkinyyQr: A.Size() == 0");
  blas_int m, lda, n, k, info, lwork;
  m = lda = A.Size(1);
  n = A.Size(2);
  if (m > n) throw LapackException("ShortLq: m > n");
  k = m;
  L.Resize(k, k);
  lwork = -1;
  double dlwork;
  gelqf(m, n, (T*) NULL, lda, (T*) NULL, &dlwork, lwork, info);
  blas_int lwork1 = (blas_int) dlwork;
  orglq(m, n, k, (T*) NULL, lda, (T*) NULL, &dlwork, lwork, info);
  blas_int lwork2 = (blas_int) dlwork;
  lwork = std::max(lwork1, lwork2);
  WorkArray<T> w(lwork + k);
  Matrix<T> work, tau;
  work.SetPtr(lwork, w.AllocWork(lwork));
  tau.SetPtr(k, w.AllocWork(k));
  gelqf(m, n, A.GetPtr(), lda, tau.GetPtr(), work.GetPtr(), lwork, info);
  T* rQ = A.GetPtr();
  T* rL = L.GetPtr();
  memset(rL, 0, k*k*sizeof(T));
  for (int c = 0; c < k; c++) {
    memcpy(rL + c, rQ + c, (m - c)*sizeof(T));
    rL += k;
    rQ += m;
  }
  orglq(m, n, k, A.GetPtr(), lda, tau.GetPtr(), work.GetPtr(), lwork, info);
}

template<typename T>
void ShortLq (const Matrix<T>& A, Matrix<T>& L, Matrix<T>& Q)
  throw (LapackException)
{
  Q = A;
  ShortLq(Q, L);
}
}
