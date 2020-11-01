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

class InverseRGreensFn : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  Matd _x;
  UInt _order;
  double _delta;

  double Eval(UInt i, UInt j) const;
};

inline double InverseRGreensFn::Eval (UInt i, UInt j) const {
  // Distance between points i and j.
  double r = 0.0;
  for (UInt k = 1; k <= 3; k++) {
    double d = _x(k, i) - _x(k, j);
    r += d*d;
  }
  r = sqrt(r + _delta);

  // If there is no softening, then I could return inf, but I prefer 0.
  if (r == 0 && _delta == 0) return 0;

  if (_order > 0)
    return 1 / pow(r, _order);
  else
    return log(r);
}

void InverseRGreensFn::Init (const KeyValueFile* kvf) throw (Exception) {
  double d;
  const Matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;

  kvf->GetDouble("delta", _delta);
  if (_delta < 0) throw Exception("delta must be >= 0.");
}

bool InverseRGreensFn::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
