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

#ifndef INCLUDE_HMMVP_COMPRESS_PRI
#define INCLUDE_HMMVP_COMPRESS_PRI

#include <vector>
#include <string>
#include "util/include/Exception.hpp"
#include "util/include/Mpi.hpp"
#include "hmmvp/include/Hmat.hpp"
#include "Hd_pri.hpp"

namespace hmmvp {
using namespace util;

// Row and column global indices for a matrix block. Indexing is 0-based.
struct MatBlock {
  UInt r0, m, c0, n;

  MatBlock();
  MatBlock(UInt r0, UInt m, UInt c0, UInt n);
  MatBlock& operator=(const hd::Block& b);
};

// Tell the user about compression of this block. Most users will likely ignore
// these data. But it's useful if you want to do element-level accuracy tricks.
struct CompressBlockInfo {
  size_t nrows, ncols;
  double tol;
  bool is_abs_tol;

  CompressBlockInfo (size_t inrows, size_t incols, bool iis_abs_tol,
                     double itol)
    : nrows(inrows), ncols(incols), tol(itol), is_abs_tol(iis_abs_tol) {}
};

// The linear algebra routines use this class to access the raw matrix.
class MatrixAccessor {
public:
  virtual ~MatrixAccessor() {}
  /* Compute B(rs,cs). Indexing starts at 0. B is preallocated.
       If rs is NULL, then blk.r0 to blk.r0 + blk.m - 1 is requested, and
     similarly for cs.
       Order B with the fastest index going down the columns.
       Return true if all is well; false if there is an error in computing the
     Green's function and you want compression to stop. */
  virtual bool Call(
    const CompressBlockInfo& cbi, const MatBlock& blk,
    const vector<UInt>* rs, const vector<UInt>* cs, double* B) const = 0;
  virtual void StartBlock (const CompressBlockInfo& cbi) const {}
  virtual void FinishBlock (const CompressBlockInfo& cbi) const {}
};

// The user sees this interface. It accounts for permutations. Indexing is
// 1-based because it's a user-side interface, and I keep all user-side
// indexing 1-based.
class GreensFn {
public:
  virtual ~GreensFn() {}
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const = 0;
  virtual void StartBlock (const CompressBlockInfo& cbi) const {}
  virtual void FinishBlock (const CompressBlockInfo& cbi) const {}
};

// Wrapper to the client's Green function to account for H-matrix permutations.
class GreensFnMa : public MatrixAccessor {
public:
  GreensFnMa(GreensFn* gf) : _gf(gf) {}

  vector<UInt>& pp () { return _pp; }
  vector<UInt>& pq () { return _pq; }
  const vector<UInt>& pp () const { return _pp; }
  const vector<UInt>& pq () const { return _pq; }

  virtual bool Call(const CompressBlockInfo& cbi, const MatBlock& blk,
                    const vector<UInt>* rs, const vector<UInt>* cs, double* B)
    const;
  virtual void StartBlock (const CompressBlockInfo& cbi) const {
    _gf->StartBlock(cbi);
  }
  virtual void FinishBlock (const CompressBlockInfo& cbi) const {
    _gf->FinishBlock(cbi);
  }

private:
  GreensFn* _gf;
  vector<UInt> _pp, _pq; // 1-based from Hd
};

// Essentially a duplicate of ACAr in Bebendorf's AHMED/ACA.h.
class AcaGfHolder;
template<typename T>
void Aca(const MatBlock& blk, MatrixAccessor& ma, double scale,
         const bool use_rel_err, const double err,
         Matrix<T>& U, Matrix<T>& V, AcaGfHolder* gh = NULL)
  throw (OutOfMemoryException, UserReqException);

// Recompression using the QR factorization.
template<typename T>
void CompressQr(Matrix<T>& U, Matrix<T>& V,
                const bool use_rel_err, const double err);

struct LraOptions {
  enum LraMethod { lra_svd, lra_aca };

  LraMethod method;
  bool is_mrem;
  double aca_tol_factor;
  UInt min_aca_size;
  bool avoid_redundant_gf_evals;
  bool call_CompressQr, allow_0rank;
  static const double qr_alpha;

  LraOptions() : method(lra_aca), aca_tol_factor(0.1), min_aca_size(16),
                 avoid_redundant_gf_evals(false), call_CompressQr(true) {}
};

template<typename T>
struct NumberBlock {
  Matrix<T> B, U, V;
};

// An LraBlock packages up everything we need for an H-matrix block. We use it
// to ferry block data without having to know its type.
class LraBlock {
public:
  virtual ~LraBlock() {}

  virtual UInt GetPrec() const = 0;
  bool IsMrem () const { return _is_mrem; }
  double GetTol () const { return _tol; }
  const MatBlock& GetBlock () const { return _b; }

  virtual void WriteToFile(FILE* fid) const = 0;
  virtual void WriteInMemory(HmatConstructor& hmc, size_t bi) const = 0;
  virtual void Pack(mpi::ByteBufferWriter& bw) const = 0;
  virtual void Unpack(mpi::ByteBufferReader& br) = 0;

protected:
  LraBlock(const MatBlock& b, bool is_mrem, double tol)
    : _b(b), _is_mrem(is_mrem), _tol(tol) {}

private:
  MatBlock _b;
  bool _is_mrem;
  double _tol;
};

template<typename T>
class TypedLraBlock : public LraBlock {
public:
  TypedLraBlock(const MatBlock& b, bool is_mrem, double tol);

  virtual UInt GetPrec() const;

  bool HaveB ()  const { return _nb.B.Size() > 0; }
  bool HaveUV () const { return _nb.U.Size() > 0; }
  UInt NcolUV () const { return _nb.U.Size(2); }

  Matrix<T>& B () { return _nb.B; }
  Matrix<T>& U () { return _nb.U; }
  Matrix<T>& V () { return _nb.V; }
  const Matrix<T>& B () const { return _nb.B; }
  const Matrix<T>& U () const { return _nb.U; }
  const Matrix<T>& V () const { return _nb.V; }

  virtual void WriteToFile(FILE* fid) const;
  virtual void WriteInMemory(HmatConstructor& hmc, size_t bi) const;
  virtual void Pack(mpi::ByteBufferWriter& bw) const;
  virtual void Unpack(mpi::ByteBufferReader& br);

private:
  NumberBlock<T> _nb;
};

// We use this factory method so that we don't have to worry about type until
// we get to the details
LraBlock* NewLraBlock(const MatBlock& b, UInt realp, bool is_mrem, double tol);

// Makes the high-level low-rank approximation decisions: whether and how to
// use an old H-matrix block, and whether to use SVD or ACA.
int ApproxByLowRank(const LraOptions& opts, LraBlock* b, MatrixAccessor& ma,
                    const LraBlock* ob = NULL)
  throw (OutOfMemoryException, UserReqException);

class HmatReader;

// Handles everything. This is exposed to the client as an opaque data type.
class Compressor {
public:
  // Matrix-wise relative error (MREM) or block-wise (BREM)?
  enum TolMethod { tm_mrem_fro, tm_mrem_abs, tm_brem_fro };

  // hd can be NULL for non-root threads.
  Compressor(const Hd* hd, GreensFn* gf) throw (Exception);
  virtual ~Compressor();

  void SetTolMethod(TolMethod tm);
  TolMethod GetTolMethod() const;
  bool IsMrem() const;
  void SetTol(double tol) throw (Exception);
  void SetBfroEstimate(double Bfro) throw (Exception);
  double GetBfroEstimate() const;
  // Use another H-matrix file to speed up compressing this one. Returns false
  // if file is incompatible.
  void UseHmatFile(const string& hmat_filename)
    throw (Exception, FileException);
  bool HaveOldHmat() const;
  double GetOldHmatBfro() const throw (Exception);
  void SetOutputLevel(UInt lev);
  UInt SetOmpNthreads(UInt nthreads);
  void AvoidRedundantGfCalls(bool flag);
  void UseCompressQr(bool flag);
  void SetPrec(UInt prec_code) throw (Exception);
  void Allow0RankBlocks(bool allow);

  double EstimateBfro()
    throw (OutOfMemoryException, UserReqException);

  void CompressToFile(const string& hmat_fn)
    throw (OutOfMemoryException, UserReqException, FileException);
  Hmat* CompressInMemory(UInt ncol, UInt max_nthreads)
    throw (OutOfMemoryException, UserReqException);

  void CompressBlockToFile(FILE* fptr, const MatBlock& b,
                           const LraBlock* ob = NULL);
  void CompressBlockInMemory(HmatConstructor& hmc, size_t bi, const MatBlock& b,
                             const LraBlock* ob = NULL);
  const string& GetFilename () const { return _filename; }
  const vector<MatBlock>& GetBlocks () const { return _blocks; }
  const vector<const MatBlock*>& GetSortedBlocks() const { return _sbs; }
  UInt GetOutputLevel() const { return _output_lvl; }
  HmatReader* GetHmatReader () { return _ohm; }
  UInt GetHmatNrows () const { return _ma.pp().size(); }
  UInt GetHmatNcols () const { return _ma.pq().size(); }

private:
  GreensFnMa _ma;

  vector<MatBlock> _blocks;
  vector<const MatBlock*> _sbs;

  string _filename; // Save to this file.

  // Use an old H-matrix to build this one.
  HmatReader* _ohm;
  bool _have_old_hmat;

  // Error tolerance data.
  UInt _prec;
  TolMethod _tm;
  double _Bfro;
  double _tol;
  bool _user_set_prec;
  bool _avoid_redundant_gf_evals, _call_CompressQr;
  bool _allow_0rank;

  UInt _output_lvl;
  UInt _omp_nthreads;

  vector<int> _ablr_stats; // Record ApproxByLowRank return states.

  static const UInt _send_sz;

private:
  void GetHdData(const Hd* hd);
  void WriteHmatHeader(FILE* fid);
  void SortBlocks();
  void CompressToFileSerial();
  void CompressToFileParallel();
  void FinalizeTol();
  void ConcatFiles();
  void RecvAndWriteFile(FILE* fid, mpi::ByteBufferReaderMpi& br, UInt pid);
  void SendFile();
  bool AmRootOrIsSerial() const;
  LraBlock* CompressBlock(const MatBlock& b, const LraBlock* ob, int& ret);
};
}

#include "Compress_inl.hpp"

#endif
