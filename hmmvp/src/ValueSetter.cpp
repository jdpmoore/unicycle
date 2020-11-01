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

#include "util/include/ValueSetter.hpp"

namespace util {
ValueSetter::ValueSetter (const mpi::ArraySegmenter* as,
                          const KeyValueFile* kvf)
  : _as(as), _kvf(kvf), _ntot(-1)
{
  _root = mpi::Root();
  _am_root = mpi::AmRoot();
}

int ValueSetter::GetNelem () const { return _n; }

void ValueSetter::InitSizes () {
  if (_ntot == -1) {
    _ntot = _as->GetNtot();
    _n = _as->GetN();
  }
}

bool ValueSetter::SetString (const string& field, string& s) {
  const string* cs;
  if (!mpi_IsTrue(_kvf->GetString(field, cs))) {
    s = string("");
    return false;
  }
  int n;
  if (_am_root) {
    s = *cs;
    n = s.size();
  }
  mpi::Bcast(&n, 1);
  _cwrk.resize(n);
  if (_am_root) memcpy(&_cwrk[0], s.c_str(), n*sizeof(char));
  mpi::Bcast(&_cwrk[0], n);
  if (!_am_root) s = string(&_cwrk[0], n);
  return true;
}

template<> bool ValueSetter::
SetArray <double>(const string& field, vector<double>& v, int factor) {
  InitSizes();
  const double* sndbuf = NULL;
  const Matd* mat;
  if (!mpi_IsTrue(_kvf->GetMatd(field, mat) &&
                  (int) mat->Size() == _ntot * factor))
    return false;
  if (_am_root) sndbuf = mat->GetPtr();
  v.resize(factor * _n);
  _as->Scatter(sndbuf, &v[0], _root, factor);
  return true;
}

template<> bool ValueSetter::
SetArray <double>(const string& field, double*& v, int factor) {
  InitSizes();
  const double* sndbuf = NULL;
  const Matd* mat;
  if (!mpi_IsTrue(_kvf->GetMatd(field, mat) &&
                  (int) mat->Size() == _ntot * factor)) {
    v = NULL;
    return false;
  }
  if (_am_root) sndbuf = mat->GetPtr();
  v = new double[factor * _n];
  _as->Scatter(sndbuf, v, _root, factor);
  return true;
}
}
