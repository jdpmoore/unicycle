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

#ifndef INCLUDE_TGF
#define INCLUDE_TGF

#include <math.h>
#include <sstream>
#include <limits>
#include "util/include/Defs.hpp"
#include "util/include/Matrix.hpp"

// Wrapper to and convenience routines for tgf, a library with code by Leslie
// Greengard, Zydrunas Gimbutas, Michael Barall, Hong Xiao. Downloaded from:
//   http://www.cims.nyu.edu/cmcl/elastostatics/tgf_2012_03_14a.zip
namespace tgf {
using namespace std;

inline void copy3 (double dst[3], const double src[3]) {
  dst[0] = src[0]; dst[1] = src[1]; dst[2] = src[2]; }
// c = a x b
inline void
calc_cross (const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0]; }
// c = b - a
inline void
calc_diff (const double a[3], const double b[3], double c[3]) {
  c[0] = b[0] - a[0];
  c[1] = b[1] - a[1];
  c[2] = b[2] - a[2]; }
inline void scale3 (double a[3], const double s) {
  a[0] *= s; a[1] *= s; a[2] *= s; }
inline double dot3 (const double a[3], const double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline double eucnorm3 (const double a[3]) { return sqrt(dot3(a, a)); }
// y += a x
inline void axpy (const double a, const double* x, double* y, const size_t n) {
  for (size_t i = 0; i < n; i++) y[i] += a*x[i]; }
inline void axpy3 (const double a, const double x[3], double y[3]) {
  y[0] += a*x[0]; y[1] += a*x[1]; y[2] += a*x[2]; }

inline void
calc_centroid (const double a[3], const double b[3], const double c[3],
               double centroid[3]) {
  static const double oth = 1.0/3;
  centroid[0] = oth*(a[0] + b[0] + c[0]);
  centroid[1] = oth*(a[1] + b[1] + c[1]);
  centroid[2] = oth*(a[2] + b[2] + c[2]);
}

// up is a unit vector that points upward along the gradient of the
// triangle. right completes the local coord system so that
//     normal = right x up.
inline void
calc_normal (const double a[3], const double b[3], const double c[3],
             double normal[3]) {
  double bma[3], cma[3];
  calc_diff(a, b, bma);
  calc_diff(a, c, cma);
  calc_cross(bma, cma, normal);
  scale3(normal, 1/eucnorm3(normal));
}
inline void calc_right (const double normal[3], double right[3]) {
  double tmp[3];
  if (max(fabs(normal[0]), fabs(normal[1])) <
      10*numeric_limits<double>::epsilon()) {
    tmp[0] = 0;
    tmp[1] = 1;
  } else {
    tmp[0] = -normal[0];
    tmp[1] = -normal[1];
  }
  tmp[2] = 0;
  calc_cross(tmp, normal, right);
  scale3(right, 1/eucnorm3(right));
}
inline void calc_up (const double normal[3], const double right[3],
                     double up[3]) {
  calc_cross(normal, right, up);
}

// normal' s v. for 3x3 symmetric tau, s = tau([1 2 3 5 6 9]).
inline double
calc_traction (const double normal[3], const double v[3], const double s[6]) {
  double ns[3];
  ns[0] = normal[0]*s[0] + normal[1]*s[1] + normal[2]*s[2];
  ns[1] = normal[0]*s[1] + normal[1]*s[3] + normal[2]*s[4];
  ns[2] = normal[0]*s[2] + normal[1]*s[4] + normal[2]*s[5];
  return dot3(ns, v);
}

inline string d3_string (const double v[3]) {
  stringstream s;
  s << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
  return s.str();
}

extern "C" void ambtgfeld_(
  const double* rlam, const double* rmu, const char* space,
  const double triangle[9], const double trinorm[3], const double obs[3],
  char* isself, const double slip[3], double u_4pi[3], double strain_4pi[9]);

// Convenience wrapper. u/strain = u/strain_4pi / (4 M_PI).
inline void eld (
  const double rlam, const double rmu,
  const double triv1[3], const double triv2[3], const double triv3[3],
  const double trinormal[3], const double obs[3], const bool isself,
  const double slip[3], double u_4pi[3], double strain_4pi[9],
  bool want_fullspace = false)
{
  char space = want_fullspace ? 'f' : 'h';
  char cisself = isself ? 't' : 'f';
  double triangle[9];
  copy3(triangle    , triv1);
  copy3(triangle + 3, triv2);
  copy3(triangle + 6, triv3);
  ambtgfeld_(&rlam, &rmu, &space, triangle, trinormal, obs, &cisself, slip,
             u_4pi, strain_4pi);
}
}

#endif
