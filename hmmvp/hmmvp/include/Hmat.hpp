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

#ifndef INCLUDE_HMMVP_HMAT
#define INCLUDE_HMMVP_HMAT

#include <vector>
#include <string>
#include "util/include/Exception.hpp"
#include "util/include/Defs.hpp"
#include "util/include/Matrix.hpp"

namespace hmmvp {
using namespace std;
using namespace util;

typedef int64 Blint;

// Initialize an H-matrix from filename. Use nthreads to compute a matrix-vector
// product. Provide space to handle MVP with arrays having up to ncol columns at
// a time.
class Hmat;
Hmat* NewHmat(const string& filename, UInt ncol = 1, UInt max_nthreads = 1,
              // Optionally read in only blocks block_idxs[0]:nb_idxs[1]-1.
              const vector<Blint>* block_idxs = 0)
  throw (FileException);

void DeleteHmat(Hmat* hm);

double HmatNormFrobenius2(const string& filename) throw (FileException);

// Perform fast matrix-vector products. This class uses just the block data
// coming from the construction of the H-matrix: in particular, the block
// cluster tree is not needed.
//   The number of OpenMP threads you use need not be in any way related to the
// number of OpenMP threads or MPI processes you used to form the H-matrix.
//   The matrix B(p,q) is approximated by A. We say that A is "permuted",
// whereas A(ip,iq), where ip = invperm(p), is "unpermuted". The unpermuted
// matrix looks like the original operator; the permuted matrix is induced by
// the cluster trees on the domain and range spaces of points. Unless you (the
// caller) want to do something very specialized, you don't need to worry about
// these permutations; leave permutation settings at their defaults.
//   Indexing begins at 0.
class Hmat {
public:
  virtual ~Hmat();

  // Number of rows.
  UInt GetM() const;
  // Number of columns.
  UInt GetN() const;
  // Number of nonzeros in the H-matrix approximation.
  UInt GetNnz() const;
  virtual UInt GetNnz(const vector<Blint>& rs, const vector<Blint>& cs) = 0;
  // sizeof(Real).
  virtual UInt GetScalarSize() const = 0;

  // Set the number of OpenMP threads to use, if OpenMP is available. Otherwise,
  // this function has (essentially) no effect. This function has to be called
  // only if you want to change the number of threads, and then only once.
  virtual UInt SetThreads(UInt nthreads) = 0;
  // Reorganize memory to get large contiguous blocks to increase cache
  // coherence. During this function call, memory use temporarily approximately
  // doubles, so this function should be called only if the memory is
  // available. If you decide to use this function, it should be called after
  // the Hmat is constructed and after each call to SetThreads. I get speedups
  // of 1 to a few tens of percents using this.
  virtual void ReorganizeMemory() = 0;

  // y = B*x, where x has ncol columns. y is allocated by the caller.
  virtual void Mvp(const float *x, float *y, UInt ncol)
    throw (Exception) = 0;
  virtual void Mvp(const double *x, double *y, UInt ncol)
    throw (Exception) = 0;
  virtual void MvpT(const float *x, float *y, UInt ncol)
    throw (Exception) = 0;
  virtual void MvpT(const double *x, double *y, UInt ncol)
    throw (Exception) = 0;
  // Subset products of the form:
  //     y = B(:,sort(cs)) x(sort(cs),:)
  //     y(sort(rs),:) = B(sort(rs),:) x
  //     y(sort(rs),:) = B(sort(rs),sort(cs)) x(sort(cs),:)
  // All of x and y must be provided, but only the subset indices are
  // meaningful. If rs = : or cs = :, then pass NULL. Passing rs = NULL, cs =
  // NULL is equivalent to calling the full Mvp.
  virtual void Mvp(const float* x, float* y, UInt ncol,
                   const vector<Blint>* rs,
                   const vector<Blint>* cs)
    throw (Exception) = 0;
  virtual void Mvp(const double* x, double* y, UInt ncol,
                   const vector<Blint>* rs,
                   const vector<Blint>* cs)
    throw (Exception) = 0;

  // Optionally disable applying permutations p and q. By default, permutations
  // are applied. If permutations are off, then (rs,cs) are relative to the
  // permuted B, not the original operator.
  void TurnOnPermute();
  void TurnOffPermute();

  // xq = Q*x
  virtual void ApplyQ(const float* x, float* xq, UInt ncol) const = 0;
  virtual void ApplyQ(const double* x, double* xq, UInt ncol) const = 0;
  // y = P'*yp
  virtual void ApplyPt(const float* yp, float* y, UInt ncol) const = 0;
  virtual void ApplyPt(const double* yp, double* y, UInt ncol) const = 0;
  virtual const Blint* GetQ() const = 0;
  virtual const Blint* GetP() const = 0;

  // Optionally save (rs,cs)-related state between calls. The caller is
  // responsible for assuring that (rs,cs) are the same between the times
  // SaveState() and ReleaseState() are called. Error checking is purposely not
  // done, as the whole idea is to speed up the MVP even more.
  //   Call SaveState() before the call to Mvp whose state should be
  // saved. Until ReleaseState() is called, subsequent calls to SaveState() have
  // no effect.
  void SaveState();
  void ReleaseState();

  // Optionally use a data structure to make subset products in which one subset
  // is small -- on the order of 10s -- faster than otherwise.
  //   Warning: The data structure can add memory overhead of ~50%, so it is not
  // always usable.
  //   I recommend against using these routines. Save/ReleaseState() already
  // speed the computation up by a lot.
  void TurnOnElementToBlockMap();
  void TurnOffElementToBlockMap();

  // These routines are not optimized. But the user should call each just once,
  // of course, since the norm does not change.
  //   Square of the Frobenius norm.
  virtual double NormFrobenius2() const = 0;
  //   1-norm.
  virtual double NormOne() const = 0;

  // Return B(rs,cs). (Unlike the Mvp routines, B(rs,cs), and not
  // B(sort(rs),sort(cs)), is returned.)
  virtual void Extract(const vector<Blint>& rs, const vector<Blint>& cs,
                       float* es) = 0;
  virtual void Extract(const vector<Blint>& rs, const vector<Blint>& cs,
                       double* es) = 0;

  // Return the (I,J,S) triple to specify the sparse matrix in Matlab
  // corresponding to this H-matrix with only the full blocks included. I and J
  // use base-1 indexing to match convention.
  //   cutoff is optional. It is a value of the sort 1e-3. Anything below
  // cutoff*max(B), where B includes only the full blocks, is discarded.
  virtual void FullBlocksIJS(
    vector<Blint>& I, vector<Blint>& J, vector<float>& S,
    float cutoff = -1.0) const = 0;
  virtual void FullBlocksIJS(
    vector<Blint>& I, vector<Blint>& J, vector<double>& S,
    double cutoff = -1.0) const = 0;

protected:
  UInt m, n; // size(B)
  UInt nnz;

  // State
  bool do_perms, save_state;
  bool state_saved;

  // Optional element-to-block map for some subset operations.
  bool do_e2b_map;

  Hmat();

private:
  friend class HmatConstructor;
};

// Build an H-matrix in memory on the fly. This is for very advanced users.
template<typename Real> class TypedHmat;
class HmatConstructor {
public:
  HmatConstructor(
    // realp = 1 for float, 2 for double.
    UInt realp,
    UInt nblocks,
    // Permutations, base-0.
    const vector<Blint>& p, const vector<Blint>& q,
    UInt ncol, UInt max_nthreads)
    throw (Exception);
  ~HmatConstructor();
  // Set block bi. r0, c0 are in H-matrix, not user, ordering (of course). The
  // caller must call the method associated with realp.
  void SetBlockB(UInt bi, UInt r0, UInt c0, const Matrix<double>& B);
  void SetBlockUV(UInt bi, UInt r0, UInt c0, const Matrix<double>& U,
                  const Matrix<double>& V);
  void SetBlockB(UInt bi, UInt r0, UInt c0, const Matrix<float>& B);
  void SetBlockUV(UInt bi, UInt r0, UInt c0, const Matrix<float>& U,
                  const Matrix<float>& V);
  // The caller must delete the H-matrix. Delete it just once, even if this
  // method is called multiple times.
  Hmat* GetHmat();
private:
  UInt _realp;
  TypedHmat<float>* _hmf;
  TypedHmat<double>* _hmd;
  bool _called_GetHmat;
};
}

#endif
