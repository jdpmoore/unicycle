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

#ifndef INCLUDE_UTIL_DEFS
#define INCLUDE_UTIL_DEFS

#include <stdlib.h>

namespace util {
  typedef size_t UInt;
}

#if __STDC_VERSION__ >= 199901L
# include <stdint.h>
// This is definitely preferred.
namespace util {
  typedef int64_t int64;
  typedef int32_t int32;
};
#else
// But I provide this for old compilers. Might run into problems with LAPACK
// integer size.
namespace util {
  typedef long long int int64;
  typedef int int32;
};
#endif

#endif
