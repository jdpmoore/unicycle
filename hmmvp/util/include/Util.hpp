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

#ifndef INCLUDE_UTIL_UTIL
#define INCLUDE_UTIL_UTIL

#include <vector>
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"

namespace util {

  template<typename T1, typename T2>
  void MatrixToVector (const Matrix<T1>& A, std::vector<T2>& v, T2 a = 0) {
    int n = A.Size();
    v.resize(n);
    const T1* pA = A.GetPtr();
    for (int i = 1; i <= n; i++) v[i] = (T2) pA[i] + a;
  }

  template<typename T1, typename T2>
  Matrix<T2>& VectorToMatrix (const std::vector<T1>& v, Matrix<T2>& A, T2 a = 0)
    throw (OutOfMemoryException)
  {
    int n = v.size();
    A.Resize(n);
    T2* pA = A.GetPtr();
    for (int i = 0; i < n; i++) pA[i] = (T2) v[i] + a;
    return A;
  }

}

#endif
