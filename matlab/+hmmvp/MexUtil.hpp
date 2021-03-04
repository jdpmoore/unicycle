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

#ifndef INCLUDE_UTILMEXUTIL
#define INCLUDE_UTILMEXUTIL

#include "mex.h"
#include "util/include/Defs.hpp"
#include "util/include/Util.hpp"

namespace util {
using namespace std;

bool GetStringv (const mxArray* ms, vector<char>& s);
bool GetString (const mxArray* ms, string& s);

template<typename T> mxArray* MatrixToMex (const Matrix<T>& m)
  throw (OutOfMemoryException)
{
  int nr = m.Size(1), nc = m.Size(2);
  mxArray* ma = mxCreateDoubleMatrix(nr, nc, mxREAL);
  if (!ma) throw OutOfMemoryException("MatrixToMex");
  double* pma = mxGetPr(ma);
  const T* pm = m.GetPtr();
  for (int i = 0, n = nr*nc; i < n; i++) pma[i] = static_cast<double>(pm[i]);
  return ma;
}

template<typename T> void MexToMatrix (const mxArray* ma, Matrix<T>& m)
  throw (OutOfMemoryException)
{
  int nr = mxGetM(ma), nc = mxGetN(ma);
  m.Resize(nr,nc);
  const double* pma = mxGetPr(ma);
  T* pm = m.GetPtr();
  for (int i = 0, n = nr*nc; i < n; i++) pm[i] = static_cast<T>(pma[i]);
}

template<typename T>
void MexToVector (const mxArray* ma, vector<T>& v, T a = 0) {
  int nr = mxGetM(ma), nc = mxGetN(ma), n = nr*nc;
  if (n == 0) return;
  v.reserve(n);
  double* pa = mxGetPr(ma);
  for (int i = 0; i < n; i++)
    v.push_back(static_cast<T>(pa[i]) + a);
}

// Could specialize MexToVector to T=double

template<typename T> mxArray* VectorToMex (const vector<T>& v, T a = 0) {
  int n = v.size();
  mxArray* ma = mxCreateDoubleMatrix(1, n, mxREAL);
  double* pma = mxGetPr(ma);
  for (int i = 0; i < n; i++) pma[i] = static_cast<double>(v[i] + a);
  return ma;
}

template<typename T>
void SetField (mxArray* s, int idx, const char* fld, const vector<T>& v,
               T idx_inc = 0) {
  mxArray* ma = VectorToMex(v, idx_inc);
  mxSetField(s, idx, fld, ma);
}

template<typename T>
void SetField (mxArray* s, int idx, const char* fld, const Matrix<T>& v) {
  mxArray* ma = MatrixToMex(v);
  mxSetField(s, idx, fld, ma);
}

template<typename T> void SetField (mxArray* s, int idx, const char* fld, T v) {
  mxArray* ma = mxCreateDoubleScalar((double)v);
  mxSetField(s, idx, fld, ma);
}

};

#endif
