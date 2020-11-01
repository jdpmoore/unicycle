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

#ifndef INCLUDE_UTIL_VALUESETTER
#define INCLUDE_UTIL_VALUESETTER

/* Given a KeyValueFile, read a field and set the part of the array for each
   processor. */

#include <vector>
#include <string>
#include "util/include/Defs.hpp"

namespace util {
using namespace std;

class KeyValueFile;
namespace mpi {
class ArraySegmenter;
}

class ValueSetter {
public:
  ValueSetter(const mpi::ArraySegmenter* as, const KeyValueFile* kvf);
  virtual ~ValueSetter() {}

  template<typename T>
  bool SetArray(const string& field, T*& v, int factor = 1);
  template<typename T>
  bool SetArray(const string& field, vector<T>& v, int factor = 1);

  template<typename T>
  void ZeroArray(T*& v, int factor = 1);

  bool SetString(const string& field, string& s);
  template<typename T>

  bool SetScalar(const string& field, T& v, T invalid = -1);

  int GetNelem() const;

private:
  const mpi::ArraySegmenter* _as;
  const KeyValueFile* _kvf;
  int _ntot, _ncomp, _n, _root;
  bool _am_root;
  vector<char> _cwrk;

  void InitSizes();
};

}

#include "ValueSetter_inl.hpp"

#endif
