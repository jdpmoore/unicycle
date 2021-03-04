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

#include "util/include/KeyValueFile.hpp"
#include "util/include/Mpi.hpp"

namespace util {

template<typename T>
bool ValueSetter::SetArray (const string& field, vector<T>& v, int factor) {
  InitSizes();
  T* sndbuf = NULL;
  const Matd* mat;
  if (!mpi_IsTrue(_kvf->GetMatd(field, mat) &&
                  (int) mat->Size() == _ntot * factor))
    return false;
  if (_am_root) {
    int fn = factor * _ntot;
    sndbuf = new T[fn];
    const double* dm = mat->GetPtr();
    for (int i = 0; i < fn; i++) sndbuf[i] = (T) dm[i];
  }
  v.resize(factor * _n);
  _as->Scatter(sndbuf, &v[0], _root, factor);
  if (_am_root) delete[] sndbuf;
  return true;
}

template<typename T>
bool ValueSetter::SetArray (const string& field, T*& v, int factor) {
  InitSizes();
  T* sndbuf = NULL;
  const Matd* mat;
  if (!mpi_IsTrue(_kvf->GetMatd(field, mat) &&
                  (int) mat->Size() == _ntot * factor)) {
    v = NULL;
    return false;
  }
  if (_am_root) {
    int fn = factor * _ntot;
    sndbuf = new T[fn];
    const double* dm = mat->GetPtr();
    for (int i = 0; i < fn; i++) sndbuf[i] = (T) dm[i];
  }
  v = new T[factor * _n];
  _as->Scatter(sndbuf, v, _root, factor);
  if (_am_root) delete[] sndbuf;
  return true;
}

template<typename T>
bool ValueSetter::SetScalar (const string& field, T& v, T invalid) {
  double d;
  if (!mpi_IsTrue(_kvf->GetDouble(field, d))) {
    v = invalid;
    return false;
  }
  if (_am_root) v = (T) d;
  mpi::Bcast(&v, 1);
  return true;
}

template<typename T>
void ValueSetter::ZeroArray (T*& v, int factor) {
  InitSizes();
  v = new T[factor * _n];
  memset(v, 0, factor * _n * sizeof(T));
}

}
