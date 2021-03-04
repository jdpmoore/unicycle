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

#ifndef INCLUDE_UTIL_LINALG
#define INCLUDE_UTIL_LINALG

#include "util/include/Defs.hpp"
#include "util/include/Matrix.hpp"
#include "util/include/WorkArray.hpp"

namespace util {
#ifdef FORTRAN_INT_4
typedef int32 blas_int;
#else
typedef int64 blas_int;
#endif

template<typename T> void gemm(
  char transa, char transb, blas_int m, blas_int n, blas_int k,
  T alpha, const T* A, blas_int lda, const T* B, blas_int ldb, T beta,
  T* C, blas_int ldc);
 
template<typename T> void gemv(
  char trans, blas_int m, blas_int n, T alpha, const T* a, blas_int lda,
  const T* x, blas_int incx, T beta, T* y, blas_int incy);

template<typename T> T dot(
  blas_int n, const T* x, blas_int incx, const T* y, blas_int incy);

template<typename T> void axpy(
  blas_int n, T a, const T* x, blas_int incx, T* y, blas_int incy);

template<typename T> void gesvd(
  char jobu, char jobvt, blas_int m, blas_int n, T* A, blas_int lda,
  T* s, T* U, blas_int ldu, T* Vt, blas_int ldvt, T* wrk, blas_int lwrk,
  blas_int& info);

template<typename T> void geqrf(
  blas_int m, blas_int n, T* A, blas_int lda, T* tau, T* work,
  blas_int lwork, blas_int& info);

template<typename T> void orgqr(
  blas_int m, blas_int n, blas_int k, T* A, blas_int lda, T* tau, T* work,
  blas_int lwork, blas_int& info);

template<typename T> void gelqf(
  blas_int m, blas_int n, T* A, blas_int lda, T* tau, T* work,
  blas_int lwork, blas_int& info);

template<typename T> void orglq(
  blas_int m, blas_int n, blas_int k, T* A, blas_int lda, T* tau, T* work,
  blas_int lwork, blas_int& info);

template<typename T> void trtrs(
  char uplo, char trans, char diag, blas_int n, blas_int nrhs, T* A,
  blas_int lda, T* B, blas_int ldb, blas_int& info);

class LapackException : public Exception {
public:
  LapackException (const std::string& msg) : Exception(msg) {}
};

template<typename T>
void Mult(const Matrix<T>& A, bool do_trans_A,
          const Matrix<T>& B, bool do_trans_B, Matrix<T>& C);

// A = R*L, with R upper tri and L lower tri.
template<typename T>
void MultRL(const Matrix<T>& R, const Matrix<T>& L, Matrix<T>& C);

// Corresponds to Matlab's [U S V] = svd(A,'econ').
template<typename T1, typename T2>
void Svd(const Matrix<T1>& A, Matrix<T2>& U, Matrix<T2>& s, Matrix<T2>& Vt)
  throw (LapackException);

// Corresponds to Matlab's [Q R] = qr(A,0) for size(A,1) >= size(A,2).
template<typename T>
void SkinnyQr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R, bool do_pr)
  throw (LapackException);
// In place version.
template<typename T>
void SkinnyQr(Matrix<T>& A, Matrix<T>& R, bool do_pr) throw (LapackException);
template<typename T>
void ShortLq(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& Q)
  throw (LapackException);
template<typename T>
void ShortLq(Matrix<T>& A, Matrix<T>& L) throw (LapackException);
}

#include "LinAlg_inl.hpp"

#endif
