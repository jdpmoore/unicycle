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

#include <algorithm>
#include "mex.h"
#include "MexUtil.hpp"
using namespace std;

namespace util {

  bool GetStringv (const mxArray* ms, vector<char>& s) {
    int strlen = mxGetNumberOfElements(ms) + 1;
    s.resize(strlen);
    if (mxGetString(ms, &s[0], strlen) != 0) return false;
    return true;
  }

  bool GetString (const mxArray* ms, string& s) {
    vector<char> vs;
    if (!GetStringv(ms, vs)) return false;
    s = string(&vs[0]);
    return true;
  }

  template<> mxArray* MatrixToMex<double> (const Matrix<double>& m)
    throw (OutOfMemoryException)
  {
    int nr = m.Size(1), nc = m.Size(2);
    mxArray* ma = mxCreateDoubleMatrix(nr, nc, mxREAL);
    if (!ma) throw OutOfMemoryException("MatrixToMex");
    memcpy(mxGetPr(ma), m.GetPtr(), nr*nc*sizeof(double));
    return ma;
  }

  template<> void MexToMatrix<double> (const mxArray* ma, Matrix<double>& m)
    throw (OutOfMemoryException)
  {
    int nr = mxGetM(ma), nc = mxGetN(ma);
    m.Resize(nr,nc);
    memcpy(m.GetPtr(), mxGetPr(ma), nr*nc*sizeof(double));
  }

};
