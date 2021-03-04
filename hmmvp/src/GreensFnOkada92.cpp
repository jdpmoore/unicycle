/* Green's function for hmmvpbuild implementing +green/stresskernel.m
   calculations from the Unicycle Matlab simulator. The basic inputs are the src
   and rcv objects.

   By the way, I (AMB) have done some work that suggests this BEM matrix may be
   very inaccurate.
     The Okada Green's function computes quantities resulting from constant slip
   over a rectangular element. If used directly for a problem on a uniform
   planar mesh, the order of accuracy for most components is 2.
     I showed (see www.stanford.edu/~ambrad/agu13.pdf for details) that if a
   planer mesh is nonuniformly discretized, then this same procedure severely
   compromises the order of accuracy, though the method still seems to be
   convergent. For the case of a planar mesh, I developed a method that restores
   the order of accuracy.
     Though I have not analyzed the case of a nonplaner mesh, my guess is that
   a nonplaner mesh, whether or not the elements are uniform, has a similar
   problem. For a mesh that has elements that are slightly overlapping or
   between which there is a slight gap, I suspect the problem is even greater.
     Therefore, the BEM matrix this file implements may be only a placeholder
   until a better method is devised. The use of an H-matrix is independent of
   these details and will be useful regardless of the details of the matrix
   entries.

   First version Feb 2014 AMB ambrad@cs.stanford.edu
*/

namespace gfb {
// Debug print in this file.
#if 0
#define bapr(fmt, ...) do {                     \
    fprintf(stderr, "GFB: " fmt, __VA_ARGS__);  \
  } while (0)
#else
#define bapr (void)
#endif

// Okada's Fortran routine DC3D.
extern "C" void dc3d_(
  const char* space, const double* ALPHA, const double* X, const double* Y,
  const double* Z, const double* DEPTH, const double* DIP, const double* AL1,
  const double* AL2, const double* AW1, const double* AW2,
  const double* DISL1, const double* DISL2, const double* DISL3,
  double* UX, double* UY, double* UZ, double* UXX, double* UYX,
  double* UZX, double* UXY, double* UYY, double* UZY, double* UXZ,
  double* UYZ, double* UZZ);

// Wrapper to it.
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

// Displacement derivatives -> stress tensor.
inline void
DuToS (const double mu, const double lambda, const double* d, double* s) {
  const double theta = d[0] + d[4] + d[8];
  s[0] = lambda*theta + 2*mu*d[0];
  s[1] = mu*(d[1] + d[3]);
  s[2] = mu*(d[2] + d[6]);
  s[3] = lambda*theta + 2*mu*d[4];
  s[4] = mu*(d[5] + d[7]);
  s[5] = lambda*theta + 2*mu*d[8];      
}

inline double dot3 (const double x[3], const double y[3]) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// Represent one element, either source or receiver.
struct Elem {
  const double *xc, *vstk, *vdip, *vnml;

  Elem (const double* xc,  const double* vstk, const double* vdip,
        const double* vnml)
    : xc(xc), vstk(vstk), vdip(vdip), vnml(vnml)
  {}
};

struct ESource : public Elem {
  double dipdeg, L, W;

  ESource (const double* xc, const double dipdeg, const double L, const double W,
           const double* vstk, const double* vdip, const double* vnml)
    : Elem(xc, vstk, vdip, vnml), dipdeg(dipdeg), L(L), W(W)
  {}

  // 1. Subtract off xc(1:2). 2. Account for strike by rotating xg, the
  // observation point in the global coord system, into xl, the obs point in the
  // element-local coord system.
  void VectorToSourceCoords (const double xg[3], double xl[3]) const {
    // v = xg - xc for x, y but not z.
    double xgc[2];
    xgc[0] = xg[0] - xc[0];
    xgc[1] = xg[1] - xc[1];
    
    // R v
    xl[0] =  vstk[0]*xgc[0] + vstk[1]*xgc[1];
    xl[1] = -vstk[1]*xgc[0] + vstk[0]*xgc[1];
    xl[2] =  xg[2];
  }

  // Undo the rotation part of VectorToSourceCoords. s is the lower triangular
  // part, stored in column-major order, of a symmetric tensor.
  void TensorToGlobalCoords (const double sl[6], double sg[6]) const {
    // R' tau
    double Rtsl[9];
    Rtsl[0] = vstk[0]*sl[0] - vstk[1]*sl[1];
    Rtsl[1] = vstk[1]*sl[0] + vstk[0]*sl[1];
    Rtsl[2] = sl[2];
    Rtsl[3] = vstk[0]*sl[1] - vstk[1]*sl[3];
    Rtsl[4] = vstk[1]*sl[1] + vstk[0]*sl[3];
    Rtsl[5] = sl[4];
    Rtsl[6] = vstk[0]*sl[2] - vstk[1]*sl[4];
    Rtsl[7] = vstk[1]*sl[2] + vstk[0]*sl[4];
    Rtsl[8] = sl[5];

    // (R' tau) R
    sg[0] = Rtsl[0]*vstk[0] - Rtsl[3]*vstk[1];
    sg[1] = Rtsl[1]*vstk[0] - Rtsl[4]*vstk[1];
    sg[2] = Rtsl[6];
    sg[3] = Rtsl[1]*vstk[1] + Rtsl[4]*vstk[0];
    sg[4] = Rtsl[2]*vstk[1] + Rtsl[5]*vstk[0];
    sg[5] = Rtsl[8];
  }
};

struct EReceiver : public Elem {
  EReceiver (const double* xc,  const double* vstk, const double* vdip,
             const double* vnml)
    : Elem(xc, vstk, vdip, vnml)
  {}

  // normal' tau v
  double CalcTraction (const double v[3], const double s[6]) const {
    double ns[3];
    ns[0] = vnml[0]*s[0] + vnml[1]*s[1] + vnml[2]*s[2];
    ns[1] = vnml[0]*s[1] + vnml[1]*s[3] + vnml[2]*s[4];
    ns[2] = vnml[0]*s[2] + vnml[1]*s[4] + vnml[2]*s[5];
    return dot3(ns, v);
  }
};

// This routine calculates the same quantities as in +green/stresskernels.m.
inline double
Dc3d (const ESource& src, const EReceiver& rcv,
      const char src_comp, const char rcv_comp, // 's' or 'd'
      const double mu, const double lambda, const double alpha) {
  double disl_strike = 0, disl_dip = 0, disl_tensile = 0;
  if (src_comp == 's') disl_strike = 1; else disl_dip = 1;

  double xr[3];
  src.VectorToSourceCoords(rcv.xc, xr);

  double u[3], du[9];
  Dc3d(alpha, xr[0], xr[1], xr[2], -src.xc[2], src.dipdeg,
       0.5*src.L, 0.5*src.L, 0.5*src.W, 0.5*src.W,
       disl_strike, disl_dip, disl_tensile, u, du);

  double sl[6];
  DuToS(mu, lambda, du, sl);
  double sg[6];
  src.TensorToGlobalCoords(sl, sg);
  
  return rcv.CalcTraction(rcv_comp == 's' ? rcv.vstk : rcv.vdip, sg);
}

// Represent the matrix to compress.
class GreensFnOkada92: public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd(double eta);
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // Subsets of the +geometry/source.m and +geometry/receiver.m class data.
  struct Patch {
    Matd xc, sv, dv, nv;
    virtual void Init(const string& kvf_fn) throw (Exception);
    void Init(const KeyValueFile* kvf) throw (Exception);
  };
  struct PSource : public Patch {
    Matd L, W, strike, dip;
    virtual void Init(const string& kvf_fn) throw (Exception);
  };
  struct PReceiver : public Patch {};

  PSource src;
  PReceiver rcv;
  double mu, lambda, alpha;
  char src_comp, rcv_comp;
  vector<ESource> esrcs;
  vector<EReceiver> ercvs;

  void ReadKvf(const KeyValueFile* kvf) throw (Exception);
};

// normal' tau v.
double CalcTraction (const double vnml[3], const double v[3], const double s[6]) {
  double ns[3];
  ns[0] = vnml[0]*s[0] + vnml[1]*s[1] + vnml[2]*s[2];
  ns[1] = vnml[0]*s[1] + vnml[1]*s[3] + vnml[2]*s[4];
  ns[2] = vnml[0]*s[2] + vnml[1]*s[4] + vnml[2]*s[5];
  return dot3(ns, v);
}

// Read the master key-value file.
void GreensFnOkada92::Init (const KeyValueFile* kvf) throw (Exception) {
  ReadKvf(kvf);

  for (size_t i = 1, N = src.L.Size(); i <= N; i++)
    esrcs.push_back(ESource(&src.xc(1,i), src.dip(i), src.L(i), src.W(i),
                            &src.sv(1,i), &src.dv(1,i), &src.nv(1,i)));
  for (size_t i = 1, N = rcv.xc.Size(2); i <= N; i++)
    ercvs.push_back(
      EReceiver(&rcv.xc(1,i), &rcv.sv(1,i), &rcv.dv(1,i), &rcv.nv(1,i)));
}

// The hierarchical decomposition uses just the element centers.
Hd* GreensFnOkada92::ComputeHd (double eta) {
  bapr("ComputeHd %ld %ld\n", rcv.xc.Size(2), src.xc.Size(2));
  return NewHd(src.xc, rcv.xc, NULL, eta);
}

// Calculate the requested matrix entries.
bool GreensFnOkada92::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Dc3d(esrcs[cs[ic]-1], ercvs[rs[ir]-1], src_comp, rcv_comp,
                  mu, lambda, alpha);
  return true;
}

// Read data from the key-value files.
#define readmat(var) do {                                               \
    if (!kvf->GetMatd(#var, m)) throw Exception("Missing " #var);       \
    Transpose(*m, var);                                                 \
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
void GreensFnOkada92::Patch::Init (const string& kvf_fn) throw (Exception) {
  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(kvf_fn))
    throw Exception("Patch::Init: Could not read " + kvf_fn);
  Init(kvf);
  DeleteKeyValueFile(kvf);
}
void GreensFnOkada92::Patch::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  readmat(xc);
  const size_t n = xc.Size(2);
  checkvarsz2(xc, 3, n);
  readmat(sv); checkvarsz2(sv, 3, n);
  readmat(dv); checkvarsz2(dv, 3, n);
  readmat(nv); checkvarsz2(nv, 3, n);  
}
void GreensFnOkada92::PSource::Init (const string& kvf_fn) throw (Exception) {
  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(kvf_fn))
    throw Exception("Source::Init: Could not read " + kvf_fn);
  Patch::Init(kvf);
  const Matd* m;
  const size_t n = xc.Size(2);
  readmat(L); checkvarsz1(L, n);
  readmat(W); checkvarsz1(W, n);
  readmat(strike); checkvarsz1(strike, n);
  readmat(dip); checkvarsz1(dip, n);
  DeleteKeyValueFile(kvf);
}
void GreensFnOkada92::ReadKvf (const KeyValueFile* kvf) throw (Exception) {
  const string* s;
  double d;
  string src_kvf_fn, rcv_kvf_fn;
  readstr(src_kvf_fn);
  readstr(rcv_kvf_fn);
  src.Init(src_kvf_fn);
  rcv.Init(rcv_kvf_fn);
  readdbl(mu);
  double nu;
  readdbl(nu);
  lambda = 2*mu*nu/(1 - 2*nu);
  alpha = (lambda + mu)/(lambda + 2*mu);
  bapr("mu %1.2e nu %1.3f lambda %1.2e alpha %1.3f\n", mu, nu, lambda, alpha);
  readchar(src_comp);
  readchar(rcv_comp);
  bapr("src_comp %c\n", src_comp);
  bapr("rcv_comp %c\n", rcv_comp);
}
#undef checkvarsz1
#undef checkvarsz2
#undef readdbl
#undef readstr
#undef readmat

#undef bapr
}
