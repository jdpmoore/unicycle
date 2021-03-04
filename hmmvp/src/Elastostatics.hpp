/* hmmvpex: Examples for hmmvp
 *   Version 0.1
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvpex is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

#ifndef INCLUDE_ELASTOSTATICS
#define INCLUDE_ELASTOSTATICS

#include <math.h>
#include <vector>
#include "util/include/Matrix.hpp"
#include "util/include/Exception.hpp"

namespace util {
using namespace std;
using namespace util;

namespace es {

inline double cosd (double a) { return cos(a*M_PI/180); }
inline double sind (double a) { return sin(a*M_PI/180); }

class LameParms {
public:
  LameParms () {}
  LameParms (double mu, double nu) { Set(mu, nu); }

  void Set (double mu,   // Shear modulus
            double nu) { // Poisson's ratio
    _mu = mu;
    _nu = nu;
    _lambda = 2*mu*nu/(1 - 2*nu);
    _alpha  = (_lambda + mu)/(_lambda + 2.0*mu);
  }

  double mu ()     const { return _mu;     }
  double nu ()     const { return _nu;     }
  double lambda () const { return _lambda; }
  double alpha ()  const { return _alpha;  }
      
private:
  double _mu, _nu, _lambda, _alpha;
};

namespace dc3 {
// Various interfaces to Y. Okada's dc3d routine.

// Declaration for the Fortran routine. The rectangle slants upward for
// positive dip in the y direction.
extern "C" void dc3d_(
  const char* space, const double* ALPHA, const double* X, const double* Y,
  const double* Z, const double* DEPTH, const double* DIP, const double* AL1,
  const double* AL2, const double* AW1, const double* AW2,
  const double* DISL1, const double* DISL2, const double* DISL3,
  double* UX, double* UY, double* UZ, double* UXX, double* UYX,
  double* UZX, double* UXY, double* UYY, double* UZY, double* UXZ,
  double* UYZ, double* UZZ);

inline void Dc3d (
  double alpha,
  double x, double y, double z, // NB: x, y are relative to the rectangle
  double depth, double dipdeg,
  double al1, double al2, double aw1, double aw2,
  double disl_strike, double disl_dip, double disl_tensile,
  double* u,  // [ux uy uz]
  double* du, // [uxx uyx uzx uxy uyy uzy uxz uyz uzz]
  bool want_fullspace = false)
{
  char wsh = 'H';
  if (want_fullspace) wsh = 'f';
  dc3d_(&wsh, &alpha,
        &x, &y, &z,
        &depth, &dipdeg, &al1, &al2, &aw1, &aw2,
        &disl_strike, &disl_dip, &disl_tensile,
        u, u+1, u+2,
        du, du+1, du+2, du+3, du+4, du+5, du+6, du+7, du+8);
}
} // namespace dc3
        
inline void DuToS (const LameParms& lp, const double* d, double* s) {
  double
    theta = d[0] + d[4] + d[8],
    mu = lp.mu(), lambda = lp.lambda();
  s[0] = lambda*theta + 2.0*mu*d[0];
  s[1] = mu*(d[1] + d[3]);
  s[2] = mu*(d[2] + d[6]);
  s[3] = lambda*theta + 2.0*mu*d[4];
  s[4] = mu*(d[5] + d[7]);
  s[5] = lambda*theta + 2.0*mu*d[8];      
}

} // namespace es
} // namespace util

#endif
