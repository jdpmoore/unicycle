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

#ifndef INCLUDE_UTIL_DEFS_PRI
#define INCLUDE_UTIL_DEFS_PRI

#include <vector>
#include "util/include/Matrix.hpp"

namespace util {

  typedef size_t UInt;
  typedef Matrix<double> Matd;
  typedef Matrix<int> Mati;
  typedef std::vector<double> Vecd;
  typedef std::vector<UInt> Vecui;

}

#if __STDC_VERSION__ >= 199901L
# include <stdint.h>
namespace util {
  typedef int64_t int64;
  typedef int32_t int32;
};
#else
namespace util {
  typedef long long int int64;
  typedef int int32;
};
#endif

#endif
