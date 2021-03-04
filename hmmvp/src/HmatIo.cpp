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

#include "util/include/Util.hpp"
#include "hmmvp/include/HmatIo.hpp"

namespace hmmvp {
  using namespace util;

  void HmatInfo (const string& filename, Blint& m, Blint& n, Blint& realp,
                 Blint& nb, double& tol)
    throw (FileException)
  {
    FILE* fid = fopen(filename.c_str(), "rb");
    if (!fid)
      throw FileException(string("HmatInfo: Can't read file ") + filename);
    ReadHmatHeader(fid, m, n, realp, nb, tol);
    fclose(fid);
  }

  void ReadHmatHeader (FILE* fid, Blint& m, Blint& n, Blint& realp, Blint& nb,
                       double& tol, Blint* rerr_method)
    throw (FileException)
  {
    FileBlint fm, fn, frealp, fnb;
    read(&fm, 1, fid);
    read(&fn, 1, fid);
    read(&frealp, 1, fid);
    fseek(fid, (fm + fn)*sizeof(FileBlint), SEEK_CUR);
    read(&fnb, 1, fid);
    m = (Blint) fm;
    n = (Blint) fn;
    realp = (Blint) frealp;
    nb = (Blint) fnb;
    FileBlint frerr_method;
    read(&frerr_method, 1, fid);
    if (rerr_method) *rerr_method = frerr_method;
    read(&tol, 1, fid);
    double blank;
    read(&blank, 1, fid);
    // This is how to construct tol for files from v0.8 to v0.13.
    if (frerr_method == 1 && blank != 14.0) tol *= blank;
  }

  void ReadHmatHeader (FILE* fid, Blint& m, Blint& n, Blint& realp, Blint& nb,
                       double& tol, vector<FileBlint>& p, vector<FileBlint>& q,
                       Blint* rerr_method)
    throw (FileException)
  {
    FileBlint fm, fn, frealp, fnb;
    read(&fm, 1, fid);
    read(&fn, 1, fid);
    read(&frealp, 1, fid);
    m = (Blint) fm;
    n = (Blint) fn;
    realp = (Blint) frealp;
    p.resize(m);
    read(&p[0], m, fid);
    q.resize(n);
    read(&q[0], n, fid);
    read(&fnb, 1, fid);
    nb = (Blint) fnb;
    FileBlint frerr_method;
    read(&frerr_method, 1, fid);
    if (rerr_method) *rerr_method = frerr_method;
    read(&tol, 1, fid);
    double blank;
    read(&blank, 1, fid);
    // This is how to construct tol for files from v0.8 to v0.13.
    if (frerr_method == 1 && blank != 14.0) tol *= blank;
  }

  void ReadHmatBlockInfo (FILE* fid, Blint realp, Blint& r0, Blint& c0,
                          Blint& m, Blint& n, Blint& rank)
    throw (FileException)
  {
    int sz = (realp == 1) ? 4 : 8;

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
      fseek(fid, (m + n)*rank*sz, SEEK_CUR);
    } else {
      rank = std::min(m, n);
      fseek(fid, m*n*sz, SEEK_CUR);
    }  
  }

}
