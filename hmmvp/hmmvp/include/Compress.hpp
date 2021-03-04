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

#ifndef INCLUDE_HMMVP_COMPRESS
#define INCLUDE_HMMVP_COMPRESS

#include <vector>
#include <string>
#include "util/include/Exception.hpp"
#include "hmmvp/include/Hd.hpp"
#include "hmmvp/include/Hmat.hpp"

namespace hmmvp {
using namespace util;
using namespace std;

// Tell the user about compression of this block. Most users will likely ignore
// these data. But it's useful if you want to do element-level accuracy tricks.
struct CompressBlockInfo {
  // Number of rows and columns in the current block.
  size_t nrows, ncols;
  // The tolerance for this block.
  double tol;
  // Whether tol is an abolute or block-wise Frobenius-norm relative error.
  bool is_abs_tol;
};

// The user implements this abstract class.
class GreensFn {
public:
  virtual ~GreensFn() {}
  // Compute B(rs,cs). Indexing starts at 1. B is preallocated. Load B with the
  // fast index on rs (like in Fortran and Matlab arrays). Return true if all is
  // well; false if there is a (catastrophic) error in computing the Green's
  // function and you want compression to stop.
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const = 0;
  // Optionally get told when Compressor is starting or finishing a block. This
  // is useful if there are block-level preprocessing steps you can do to speed
  // up the computation.
  virtual void StartBlock (const CompressBlockInfo& cbi) const {}
  virtual void FinishBlock (const CompressBlockInfo& cbi) const {}
};

// Build an H-matrix.
class Compressor {
private:
  virtual ~Compressor();

public:

  // Error tolerance.

  enum TolMethod {
    // MREM with error bound tol ||B||_F. This is the recommended method but is
    // not the default because ||B||_F must be set. It builds the H-matrix
    // approximation A to B so that ||B - A||_F <= tol ||B||_F.
    tm_mrem_fro,
    // MREM with error bound tol. This is the default method. It builds the
    // H-matrix so that ||B - A||_F <= tol.
    tm_mrem_abs,
    // BREM with error bound tol ||B||_F. The bound ||B - A||_F <= tol ||B||_F
    // is the same as with tm_mrem_fro, but it also enforces a block-level
    // relative error from which some applications can benefit. It usually is
    // costlier, by factors of 1.5 to 4, to use tm_brem_fro than tm_mrem_fro.
    tm_brem_fro
  };
  void SetTolMethod(TolMethod tm);
  TolMethod GetTolMethod() const;
  // If tm_mrem_fro is used, set ||B||_F. An approximation within a few 10s of
  // percent is enough in practice.
  void SetBfroEstimate(double Bfro) throw (Exception);
  // One good way to get an estimate of ||B||_F is to call this method. It only
  // works well when the matrix A is dominant on and near its diagonal, and it
  // may work terribly in other cases.
  double EstimateBfro()
    throw (OutOfMemoryException, UserReqException);
  // See what value for ||B||_F is being used.
  double GetBfroEstimate() const;
  // Using one of the tm_mrem_* methods?
  bool IsMrem() const;
  // Set tol. Default is 1e-6.
  void SetTol(double tol) throw (Exception);

  // Using an old H-matrix.

  // Use another H-matrix file to speed up compressing this one. Returns false
  // if file is incompatible.
  void UseHmatFile(const string& hmat_filename)
    throw (Exception, FileException);
  // Are we using a compatible old H-matrix to speed up compression?
  bool HaveOldHmat() const;
  // An even better way to get an estimate of ||B||_F is to calculate it from an
  // old H-matrix. Evan an inaccurate H-matrix will yield a very good estimate
  // for ||B||_F, and dominance near the diagonal is unnecessary.
  double GetOldHmatBfro() const throw (Exception);

  // Tweaks.

  // 0 for little output, > 0 for more.
  void SetOutputLevel(UInt lev);
  // OpenMP parallelization is permitted only if a single Compressor exists at a
  // time. Otherwise, nthreads must be 1. If using MPI, this method does
  // nothing.
  UInt SetOmpNthreads(UInt nthreads);
  // Reduce the number of Green's function calls at the expense of some memory
  // overhead.
  void AvoidRedundantGfCalls(bool flag);
  // There's almost never a reason to call this method.
  void UseCompressQr(bool flag);
  
  // Run.

  // Save the H-matrix to a file.
  void CompressToFile(const string& hmat_fn)
    throw (OutOfMemoryException, UserReqException, FileException);
  // Construct the H-matrix in memory. This is for advanced users.
  Hmat* CompressInMemory(UInt ncol, UInt max_nthreads)
    throw (OutOfMemoryException, UserReqException);
  
  // This shouldn't actually be called unless I'm doing a convergence test.
  void SetPrec(UInt prec_code) throw (Exception);
  // This shouldn't be called either except for a tol test. In practice I
  // recommend using a rank-1 block even when a rank-0 (all 0) block meets the
  // allocated tolerance.
  void Allow0RankBlocks(bool allow);

private:
  Compressor();
  Compressor(const Compressor&);
  Compressor& operator=(const Compressor&);
};

// Construct and delete a Compressor.
Compressor* NewCompressor(const Hd* hd, GreensFn* gf) throw (Exception);
void DeleteCompressor(Compressor* c);
}

#endif
