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

#ifndef INCLUDE_HMMVP_HMATIO
#define INCLUDE_HMMVP_HMATIO

#include "util/include/Defs.hpp"
#include "hmmvp/include/Hmat.hpp"

// Interface to read an H-matrix file. One good use for these routines is to
// implement a matrix-vector product specialized to your computer setup.

namespace hmmvp {

typedef int64 FileBlint;

// Code for single (1) and double precision (2) storage.
template<typename T>
unsigned char GetPrecCode () { return sizeof(T) == sizeof(float) ? 1 : 2; }

// Given a file pointer from fopen, read the size of the H-matrix (m by n), the
// precision of the storage (realp), the number of blocks (nb), the tolerance
// (tol: multiple interpretations, depending on error method), base-0 row (p)
// and column permutations (q) that permute the original matrix into
// block-low-rank revealing form, and optionally the error method (1 for MREM, 0
// for BREM). Generally, the caller does not care about tol and rerr_method.
void ReadHmatHeader(FILE* fid, Blint& m, Blint& n, Blint& realp, Blint& nb,
                    double& tol, vector<FileBlint>& p, vector<FileBlint>& q,
                    Blint* rerr_method = NULL)
  throw (FileException);

// Disregard p and q.
void ReadHmatHeader(FILE* fid, Blint& m, Blint& n, Blint& realp, Blint& nb,
                    double& tol, Blint* rerr_method = NULL)
  throw (FileException);

// After reading the header by either of the previous functions, sequentially
// read each block. All data refer to the permuted matrix B, not the
// original. r0 and c0 are the base-0 row and column start of the block. m
// (rows) and n (columns) are the block size. rank is the block's rank. On
// entrance, B, U, and Vt should be uninitialized pointers. On exit, either of
// these situations holds: (1) B is NULL and U and Vt contain the block low-rank
// approximation U*Vt, or (2) U and Vt are NULL and B contains the exact block
// data.
template<typename T>
void ReadHmatBlock(FILE* fid, Blint& r0, Blint& c0, Blint& m, Blint& n,
                   Blint& rank, T*& B, T*& U, T*& Vt)
  throw (FileException);

// Disregard the data.
void ReadHmatBlockInfo(FILE* fid, Blint realp, Blint& r0, Blint& c0, Blint& m,
                       Blint& n, Blint& rank)
  throw (FileException);

// Also disregard the precision.
template<typename Real>
void ReadHmatBlockInfo(FILE* fid, Blint& r0, Blint& c0, Blint& m, Blint& n,
                       Blint& rank)
  throw (FileException);

// Convenience routine. Return the size and precision of the H-matrix stored in
// the file filename.
void HmatInfo(const string& filename, Blint& m, Blint& n, Blint& realp,
              Blint& nb, double& tol)
  throw (FileException);

}

#include "HmatIo_inl.hpp"

#endif
