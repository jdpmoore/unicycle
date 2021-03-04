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

#ifndef INCLUDE_HMMVP_HD_PRI
#define INCLUDE_HMMVP_HD_PRI

#include <assert.h>
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"
#include "util/include/Defs.hpp"

// Create an H-matrix from B. P'BQ has contiguous blocks.

namespace hmmvp {
using namespace std;
using namespace util;

typedef vector<UInt> Vecui;
typedef Matrix<double> Matd;
  
// Row and column global indices for a matrix block. Indexing is 1-based.
namespace hd {
struct Block {
  UInt r0, m, c0, n;

  Block(UInt r0, UInt m, UInt c0, UInt n);
  Block(const Vecui& p, const Vecui& q, const Vecui& r, const Vecui& c);
};
}

// Block cluster tree, represented by just the node. struct Block is a
// leaf. Later, we make a list<Block>, and block_idx is an index into this list
// if BNode is a leaf, else -1. struct Block is actually redundant relative to a
// full block cluster tree, but for a matrix-vector product, the tree is not
// necessary.
struct BNode {
  BNode* chld[4];
  int block_idx;
  UInt r0, m, c0, n;

  BNode(UInt r0, UInt m, UInt c0, UInt n);
  ~BNode();
  UInt Nchldrn() const;
  UInt Nnodes() const;
};

class Hd {
public:
  static const UInt cluster_tree_min_points;

  Hd() : nblocks(0) {}
  virtual ~Hd();

  // Permutations use base-1 indexing.
  virtual void Permutations(Vecui& p, Vecui& q) const = 0;

  // Methods to access matrix blocks.
  typedef vector<hd::Block>::const_iterator iterator;
  iterator Begin() const;
  iterator End() const;
  vector<hd::Block>::size_type NbrBlocks() const;

  // Handle periodic boundary conditions. dim is a base-1 idx.
  void SetPeriodicBoundary(int dim, double lo, double hi);
  // Distance squared, accounting for possible periodic BC.
  double Distance2(const double* p1, const double* p2, int nd) const;

protected:
  // Matrix blocks
  vector<hd::Block> blocks;
  UInt nblocks;
  // Periodic boundaries
  Matd _pbounds;
  vector<bool> _is_periodic;
  vector<double> _pwidths;
};

class Node;
class HdCcCD;

// H-matrix hierarchical decomposition, cell-centered.
class HdCc : public Hd {
public:
  HdCc();
  virtual ~HdCc() {}

  virtual void Compute() = 0;

  void Opt_eta(double eta) throw (Exception);
  // Cluster splitting method. Default is to split a method along the major axis
  // of the cluster. Alternatively, split along axis directions only.
  enum SplitMethod { splitMajorAxis, splitAxisAligned };
  void Opt_split_method(SplitMethod sm) { split_method = sm; }
  // Method to determine whether two clusters are far enough apart. Default is
  // distCentroid.
  enum DistanceMethod { distCentroid, distPointwise, distAxisAligned };
  void Opt_distance_method(DistanceMethod sm) { dist_method = sm; }
    
protected:
  bool Admissible(const HdCcCD& n1, const HdCcCD& n2) const;

  // eta as used in the standard admissibility condition.
  double eta;
  SplitMethod split_method;
  DistanceMethod dist_method;
};  

class Node;

// Use this class if the domain (also called influence or source) points are the
// same as the range points. These points are in the 3xN matrix A.
class HdCcSym : public HdCc {
public:
  HdCcSym(const Matd& A, const Matd* wgt = 0);
  virtual void Compute();
  BNode* Compute(bool want_block_cluster_tree);

  // Permutation matrix so that the blocks are contiguous. Indexing is
  // 1-based. P'BP has contiguous blocks.
  virtual void Permutations(Vecui& p, Vecui& q) const
  { p = permp; q = permp; }

private:
  void Compute(const vector<const Node*>& ns1,
               const vector<const Node*>& ns2);
  BNode* Compute(const Node* n1, const Node* n2);

private:
  // Domain and range meshs
  Matd A;
  Vecui permp, permpi;
};

// Use this class if the domain (D) and range (R) points are different or one is
// a subset of the other. This class returns the same results as HdCcSym (with
// permutations p == q) if D == R; but this class takes ~25% longer in this
// case.
class HdCcNonsym : public HdCc {
public:
  HdCcNonsym(const Matd& domain, const Matd& range, const Matd* wgt = 0);
  virtual void Compute();

  // Permutation matrices so that the blocks are contiguous. Indexing is
  // 1-based. P'BQ has contiguous blocks.
  virtual void Permutations(Vecui& p, Vecui& q) const
  { p = permp; q = permq; }

private:
  void Compute(const vector<const Node*>& ns1,
               const vector<const Node*>& ns2);

private:
  // Domain and range meshs
  Matd D, R;
  Vecui permp, permpi, permq, permqi;
};

// For use in reading a spatial decomposition from a file.
class HdData : public Hd {
public:
  virtual void Permutations(Vecui& p, Vecui& q) const
  { p = permp; q = permq; }

  Vecui& Permp() { return permp; }
  Vecui& Permq() { return permq; }
  void AppendBlock(const hd::Block& b) { blocks.push_back(b); nblocks++; }

private:
  Vecui permp, permq;
};

}

#endif
