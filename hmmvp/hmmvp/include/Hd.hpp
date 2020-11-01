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

#ifndef INCLUDE_HMMVP_HD
#define INCLUDE_HMMVP_HD

#include <vector>
#include "util/include/Defs.hpp"
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"

namespace hmmvp {
using namespace util;

// Row and column global indices for a matrix block. Indexing is 1-based.
namespace hd {
struct Block { UInt r0, m, c0, n; };
}

class Hd {
private:
  virtual ~Hd();

public:
  // Permutations use base-1 indexing.
  virtual void Permutations(std::vector<UInt>& p, std::vector<UInt>& q)
    const = 0;

  // Methods to access matrix blocks.
  typedef std::vector<hd::Block>::const_iterator iterator;
  iterator Begin() const;
  iterator End() const;
  std::vector<hd::Block>::size_type NbrBlocks() const;

private:
  Hd();
  Hd(const Hd&);
  Hd& operator=(const Hd&);
};

/* R and D are 3xN matrices. D corresponds to the influence (source) points;
   these are the columns of the matrix. R corresponds to the influenced points,
   which are the rows. In many cases, D == R, in which case you do not have to
   specify R. R can be empty or not provided if R == D.
     Optionally provide a 2x3 matrix setting the coordinates of axis-aligned
   periodic boundaries. For example, if the domain is periodic in the x
   direction, and the primary domain extends from x = 2 to x = 3, then set
   periodic_boundaries(1, 1) = 2 and periodic_boundaries(1, 2) = 3. Set
   periodic_boundaries(1, i) >= periodic_boundaries(2, i) if dimension i
   (base-1) is not periodic. */
Hd* NewHd(const Matrix<double>& D,
          const Matrix<double>* periodic_boundaries = NULL,
          double eta = 3);
Hd* NewHd(const Matrix<double>& D, const Matrix<double>& R,
          const Matrix<double>* periodic_boundaries = NULL,
          double eta = 3);

void WriteHd(const Hd* hd, const std::string& hd_filename)
  throw (FileException);

Hd* NewHd(const std::string& hd_filename)
  throw (Exception, FileException);

void DeleteHd(Hd* hd);

#ifdef TESTING_AND_ANALYSIS
/* Splitting along the major axis is the general method and gives perfectly fine
   results even when splitting along the major axes is applicable. But for
   precise tol analysis, I want to split along major axes, so I use these: */
Hd* NewHdAxisAligned(const Matrix<double>& D, const Matrix<double>* pb,
                     double eta = 3);

#endif
}

#endif
