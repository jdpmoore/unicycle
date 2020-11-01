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

/* Example Green's function implementation that wraps tgf. Downloaded from
       http://www.cims.nyu.edu/cmcl/elastostatics/tgf_2012_03_14a.zip
   on 17 Feb 2014. Some of the methods described in:
       Gimbutas, Zydrunas, et al. "On the Calculation of Displacement, Stress,
       and Strain Induced by Triangular Dislocations." Bulletin of the
       Seismological Society of America 102.6 (2012): 2776-2780.

   Disclaimer: This file is meant to be an example of how to insert new
   ImplGreensFns into hmmvp. It is not very well tested.

   'hmmvpbuild help compress' gives the following:

    tgf: Wrapper to tgf [Gibmutas et al], which is used to compute traction due
         to constant displacement on a triangle in a homogeneous elastic
         fullspace or halfspace.
      vertices [3xNv array]: Triangle vertices. Nv is the number of vertices.
      triangles [3xNt array]: Index (starting at 1) into vertices describing the
        triangles.
      mu [scalar]: Shear modulus.
      nu [scalar]: Poisson's ratio.
      disl [3-vector]: Dislocation vector in local coordinates. x component
        is strike-slip, y is dip-slip, z is tensile.
      rcv [3-vector]: Receiver traction vector in local coordinates.
      want_fullspace [scalar, optional, default 0]: Non-0 if a fullspace, rather
        than a halfspace, is desired.
*/

#include "Elastostatics.hpp"
#include "tgf.hpp"

namespace tgf {
using namespace util;
using namespace es;
using namespace dc3;

typedef Matrix<double> Matd;

class GreensFn : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd(double eta);
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  void calc_tri_data(const double disl_lcl[3], const double rcv_lcl[3]);

  Matd v, normal, centroid, src, rcv;
  Matrix<int> tri;
  LameParms lp;
  bool want_fullspace;
};

void GreensFn::
calc_tri_data (const double disl_lcl[3], const double rcv_lcl[3]) {
  const size_t n = tri.Size(2);
  normal.Reshape(3, n);
  centroid.Reshape(3, n);
  src.Reshape(3, n);
  rcv.Reshape(3, n);
  for (size_t i = 1; i <= n; i++) {
    calc_centroid(&v(1,tri(1,i)), &v(1,tri(2,i)), &v(1,tri(3,i)),
                  &centroid(1,i));

    // Local coord system in global coord.
    double along_strike[3], along_dip[3];
    calc_normal(&v(1,tri(1,i)), &v(1,tri(2,i)), &v(1,tri(3,i)), &normal(1,i));
    calc_right(&normal(1,i), along_strike);
    calc_up(&normal(1,i), along_strike, along_dip);

    // Local dislocation in global coord.
    axpy3(disl_lcl[0], along_strike, &src(1,i));
    axpy3(disl_lcl[1], along_dip, &src(1,i));
    axpy3(disl_lcl[2], &normal(1,i), &src(1,i));

    // Local receiver vector in global coord.
    axpy3(rcv_lcl[0], along_strike, &rcv(1,i));
    axpy3(rcv_lcl[1], along_dip, &rcv(1,i));
    axpy3(rcv_lcl[2], &normal(1,i), &rcv(1,i));
  }
}

Hd* GreensFn::ComputeHd (double eta) {
  return NewHd(centroid, centroid, NULL, eta);
}

bool GreensFn::Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
                     const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++) {
    const UInt c = cs[ic];
    const int* tri_ic = &tri(1,c);
    for (UInt ir = 0; ir < rs.size(); ir++, k++) {
      const UInt r = rs[ir];
      double u[3], strain[9], stress[6];
      eld(lp.lambda(), lp.mu(),
          &v(1,tri_ic[0]), &v(1,tri_ic[1]), &v(1,tri_ic[2]),
          &normal(1,c), &centroid(1,r), c == r, &src(1,c), u, strain,
          want_fullspace);
      DuToS(lp, strain, stress);
      B[k] = calc_traction(&normal(1,r), &rcv(1,r), stress)/(4*M_PI);
      assert(!isnan(B[k]));
    }
  }
  return true;
}

// Read data from the key-value files.
#define readmat(name, var) do {                                         \
    if (!kvf->GetMatd(#name, m)) throw Exception("Missing " #name);     \
    var = *m;                                                           \
  } while (0)
#define readstr(var) do {                                               \
    if (!kvf->GetString(#var, s)) throw Exception("Missing " #var);     \
    var = *s;                                                           \
  } while (0)
#define readchar(var) do {                                              \
    if (!kvf->GetString(#var, s)) throw Exception("Missing " #var);     \
    if (s->size() != 1) throw Exception("Must be a char: " #var);       \
    var = (*s)[0];                                                      \
  } while (0)
#define readdbl(var) do {                                               \
    if (!kvf->GetDouble(#var, d)) throw Exception("Missing " #var);     \
    var = d;                                                            \
  } while (0)
#define checkvarsz2(var, m, n) do {             \
    if (var.Size(1) != m || var.Size(2) != n)   \
      throw Exception("Wrong size: " #var);     \
  } while (0)
#define checkvarsz1(var, m) do {                \
    if (var.Size() != m)                        \
      throw Exception("Wrong size: " #var);     \
  } while (0)

void GreensFn::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  double d;

  double mu, nu;
  readdbl(mu);
  readdbl(nu);
  lp.Set(mu, nu);

  readmat(vertices, v);
  const size_t nv = v.Size(2);
  checkvarsz2(v, 3, nv);

  Matrix<double> dtri;
  readmat(triangles, dtri);
  const size_t ntri = dtri.Size(2);
  checkvarsz2(dtri, 3, ntri);
  tri.Resize(3, ntri);
  for (size_t j = 1; j <= ntri; j++)
    for (size_t i = 1; i <= 3; i++) {
      tri(i,j) = (int) dtri(i,j);
      if (tri(i,j) < 1 || tri(i,j) > (int) nv)
        throw Exception("tri entry is out of range.");
    }

  Matrix<double> disl(3), rcv(3);
  readmat(disl, disl); checkvarsz1(disl, 3);
  readmat(rcv, rcv); checkvarsz1(rcv, 3);

  want_fullspace = false;
  if (kvf->GetDouble("want_fullspace", d)) want_fullspace = (bool) d;

  calc_tri_data(disl.GetPtr(), rcv.GetPtr());
}

#undef checkvarsz1
#undef checkvarsz2
#undef readdbl
#undef readstr
#undef readmat

#undef bapr
}
