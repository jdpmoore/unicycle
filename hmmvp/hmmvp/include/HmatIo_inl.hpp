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

#ifndef INCLUDE_HMMVP_HMATIO_INL
#define INCLUDE_HMMVP_HMATIO_INL

#include "util/include/Util.hpp"

namespace hmmvp {
  using namespace util;
  
  template<typename Real>
  void ReadHmatBlockInfo (FILE* fid, Blint& r0, Blint& c0, Blint& m, Blint& n,
                          Blint& rank)
    throw (FileException)
  {
    Blint realp = sizeof(Real) == 4 ? 1 : 2;
    ReadHmatBlockInfo(fid, realp, r0, c0, m, n, rank);
  }

  template<typename Real>
  void ReadHmatBlock (FILE* fid, Blint& r0, Blint& c0, Blint& m, Blint& n,
                      Blint& rank, Real*& B, Real*& U, Real*& Vt)
    throw (FileException)
  {
    FileBlint fr0, fc0, fm, fn;
    read(&fr0, 1, fid);
    read(&fm, 1, fid);
    read(&fc0, 1, fid);
    read(&fn, 1, fid);
    r0 = (Blint) fr0;
    c0 = (Blint) fc0;
    m = (Blint) fm;
    n = (Blint) fn;

    char code;
    read(&code, 1, fid);

    if (code == 'U') {
      FileBlint frank;
      read(&frank, 1, fid);
      rank = (Blint) frank;
      Vt = new Real[n*rank];
      U = new Real[m*rank];
      read(U, m*rank, fid);
      read(Vt, n*rank, fid);
      B = 0;
    } else {
      rank = std::min(m, n);
      B = new Real[m*n];
      read(B, m*n, fid);
      U = 0;
      Vt = 0;
    }
  }

}

#endif
