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

#ifndef INCLUDE_UTIL_IO
#define INCLUDE_UTIL_IO

#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include "util/include/Exception.hpp"

namespace util {

  template<typename T>
  void read (T* v, size_t sz, FILE* fid) throw (FileException) {
    size_t nread;
    if (!fid) throw FileException("read: fid is null");
    if ((nread = fread(v, sizeof(T), sz, fid)) != sz) {
      std::stringstream s;
      s << "read: nread = " << nread << " sz = " << sz;
      throw FileException(s.str());
    }
  }

  template<typename T>
  void write (const T* v, size_t sz, FILE* fid) throw (FileException) {
    size_t nwrite;
    if (!fid) throw FileException("write: fid is null");
    if ((nwrite = fwrite(v, sizeof(T), sz, fid)) != sz) {
      std::stringstream s;
      s << "write: nwrite = " << nwrite << " sz = " << sz;
      throw FileException(s.str());
    }
  }

  template<typename T>
  void Write (const std::vector<T>& v, FILE* fid) throw (FileException) {
    size_t nv;
    nv = v.size();
    write(&nv, 1, fid);
    write(&v[0], nv, fid);
  }

  template<typename T> void Read (std::vector<T>& v, FILE* fid) {
    size_t nv;
    read(&nv, 1, fid);
    v.resize(nv);
    read(&v[0], nv, fid);
  }

}

#endif
