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

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include "util/include/Util.hpp"
#include "util/include/LinAlg.hpp"
#include "Hd_pri.hpp"

namespace hmmvp {
const UInt Hd::cluster_tree_min_points = 3;

// Inverse of a permutation vector p. delta = 1 is for 1-based indexing.
inline void InvPermutation (const Vecui& p, Vecui& ip, UInt delta = 1) {
  UInt n = p.size();
  ip.resize(n);
  for (UInt i = 0; i < n; i++) ip[p[i] - delta] = i + delta;
}

// Find the minimum index in the block in the permuted matrix.
UInt MinBlockIdx (const Vecui& ip, const Vecui& is) {
  UInt n = is.size(), mi = ip.size();
  for (UInt i = 0; i < n; i++) {
    UInt ipi = ip[is[i]-1];
    if (ipi < mi) mi = ipi;
  }
  return mi;
}

hd::Block::Block (UInt r0_, UInt m_, UInt c0_, UInt n_)
  : r0(r0_), m(m_), c0(c0_), n(n_)
{
  assert(r0 >= 1 && c0 >= 1 && m >= 1 && n >= 1);
}

// p is the inverse permuation matrix; r and c are indices into the unpermuted
// matrix.
hd::Block::
Block (const Vecui& ip, const Vecui& iq, const Vecui& r, const Vecui& c) {
  r0 = MinBlockIdx(ip, r);
  m = r.size();
  c0 = MinBlockIdx(iq, c);
  n = c.size();
}

BNode::BNode(UInt r0, UInt m, UInt c0, UInt n)
  : block_idx(-1), r0(r0), m(m), c0(c0), n(n)
{
  for (int i = 0; i < 4; i++) chld[i] = 0;
}

BNode::~BNode () {
  for (int i = 0; i < 4; i++)
    if (chld[i]) delete chld[i];
}

UInt BNode::Nnodes () const {
  UInt n = 1;
  for (int i = 0; i < 4; i++)
    if (chld[i]) n += chld[i]->Nnodes();
  return n;
}

UInt BNode::Nchldrn () const {
  UInt nc = 0;
  for (int i = 0; i < 4; i++)
    if (chld[i]) nc++;
  return nc;
}

// B = A(:,cs)
void MatrixSubCols (const Matd& A, const Vecui& cs, Matd& B) {
  int m = A.Size(1), Bn = cs.size();
  B.Resize(m, Bn);
  int Aos = 0, Bos = 0;
  for (int j = 0; j < Bn; j++) {
    Aos = m*(cs[j] - 1);
    memcpy(B.GetPtr() + Bos, A.GetPtr() + Aos, m*sizeof(double));
    Bos += m;
  }
}

void MajorAxis (const Matd& A, Matd& u) {
  // Centroid of cluster
  int m = A.Size(1), n = A.Size(2);
  Matd ctr(m);
  ctr.Bzero();
  for (int j = 1; j <= n; j++)
    for (int i = 1; i <= m; i++)
      ctr(i) += A(i,j);
  ctr *= 1.0 / n;
  // Subtract off centroid
  Matd B = A;
  for (int j = 1; j <= n; j++)
    for (int i = 1; i <= m; i++)
      B(i,j) -= ctr(i);

  Matd U, s, Vt;
  Svd(B, U, s, Vt);
  u.Resize(m);
  memcpy(u.GetPtr(), U.GetPtr(), u.Size(1)*sizeof(double));
}

template<typename T> ostream& operator<< (ostream& os, const vector<T>& v) {
  os << "[ ";
  for (unsigned int i = 0; i < v.size(); i++) os << v[i] << " ";
  os << "]";
  return os;
}

// ---- ClusterTree

// Abtract class for holding the data in each cluster. The data are used to
// determine whether two clusters can form a block (are "admissible").
class ClusterData {
public:
  virtual ~ClusterData() {};
  virtual void Init(const Vecui& idxs, const Matd& A, const Matd* wgt) = 0;
};

class ClusterDataFactory {
public:
  virtual ~ClusterDataFactory () {}
  virtual ClusterData* CreateClusterData() = 0;
};

// Node in the cluster tree.
class Node {
public:
  Node () : cd(0) {}
  ~Node();
  
  // Node IDs are set so that for nodes n1 and n2 at level L1, if n1.Id() <
  // n2.Id(), then for every descendant d1 of n1 and d2 of n2, d1.Id() <
  // d2.Id().
  int Id () const { return id; }
  const Vecui& Idxs () const { return idxs; }
  ClusterData& CD () const { return *cd; }
  const vector<const Node*>& Chldrn () const { return const_chldrn; }

private:
  friend class ClusterTree;

  Vecui idxs;
  int id;
  ClusterData* cd;
  // TODO: We use just two children. Should make this Node* chldrn[2].
  vector<Node*> chldrn;
  vector<const Node*> const_chldrn;
};

// Cluster tree for general n-dimensional nodes.
class ClusterTree {
public:
  ClusterTree(ClusterDataFactory* cdf, HdCc::SplitMethod sm, const Matd* A,
              const Matd* wgt = 0);
  virtual ~ClusterTree();

  void Opts_min_points (int set) { min_points = set; }

  void Compute(bool setid = true);

  Node& Root () const { return *root; }
  // Leaves are ordered according to a depth-first search.
  void Leaves(vector<const Node*>& lvs) const;

  // Get the permutation of the indices of a symmetric matrix implied by this
  // cluster tree.
  void Permutation(Vecui& p) const;

protected:
  virtual void Compute(Node* node);

private:
  int SetId(Node* node, int id);
  void Leaves(vector<const Node*>& lvs, Node* node) const;
  
protected:
  ClusterDataFactory* cdf;
  // Split on major axis or aligned with axes
  HdCc::SplitMethod split_method;
  // nd x N array of points
  const Matd& A;
  Node* root;
  // Minimum number of points in a cluster
  int min_points;
};

Node::~Node () {
  if (cd) delete cd;
  for (unsigned int i = 0; i < chldrn.size(); i++)
    delete chldrn[i];
}

ClusterTree::
ClusterTree (ClusterDataFactory* cdf, HdCc::SplitMethod sm, const Matd* A,
             const Matd* wgt)
  : cdf(cdf), split_method(sm), A(*A), root(0),
    min_points(Hd::cluster_tree_min_points)
{}

ClusterTree::~ClusterTree () {
  if (root) delete root;
}

void ClusterTree::Compute (bool setid) {
  // Root node
  int n = A.Size(2);
  if (n == 0) return;
  root = new Node();
  Vecui& is = root->idxs;
  is.resize(n);
  for (int i = 1; i <= n; i++) is[i-1] = i;
  ClusterData* cd = cdf->CreateClusterData();
  cd->Init(is, A, 0);
  root->cd = cd;

  // Build tree
  Compute(root);

  if (setid) {
    // Give each node an id
    SetId(root, 0);
  }
}

void ClusterTree::Compute (Node* node) {
  Vecui& is = node->idxs;
  int n = is.size();
  if (n <= min_points) return;

  // New nodes
  for (int i = 0; i < 2; i++) {
    node->chldrn.push_back(new Node());
    node->const_chldrn.push_back(node->chldrn[i]);
  }
  Vecui &is1 = node->chldrn[0]->idxs, &is2 = node->chldrn[1]->idxs;
  is1.reserve(n);
  is2.reserve(n);

  if (split_method == HdCc::splitMajorAxis) {
    Matd Ais, u, p(n);
    // Get major axis of A(:,is)
    MatrixSubCols(A, is, Ais);
    MajorAxis(Ais, u);
    // Project onto major axis: p = u'*A
    u.Reshape(1, u.Size());
    mmp(u, Ais, p);
    // Split based on projection
    double mip = p(1), map = p(1);
    for (int i = 2; i <= n; i++) {
      double pi = p(i);
      if (pi < mip) mip = pi;
      else if (pi > map) map = pi;
    }
    double dp = map - mip, cut = mip + dp/2.0;
    for (int i = 1; i <= n; i++) {
      double pi = p(i);
      if (pi <= cut) is1.push_back(is[i-1]);
      else is2.push_back(is[i-1]);
    }
  } else { // HdCc::splitAxisAligned
    int nd = A.Size(1);
    // Get axis-aligned box limits
    Matd szs(nd,2);
    for (Vecui::size_type i = 0; i < is.size(); i++) {
      int k = is[i];
      for (int j = 1; j <= nd; j++)
        if (i == 0) {
          szs(j,1) = szs(j,2) = A(j,k);
        } else {
          if (A(j,k) < szs(j,1)) szs(j,1) = A(j,k);
          if (A(j,k) > szs(j,2)) szs(j,2) = A(j,k);
        }
    }
    // Find dimension having maximum box extent
    int kd = 1;
    double maxlen = szs(1,2) - szs(1,1);
    for (int i = 2; i <= nd; i++) {
      double d = szs(i,2) - szs(i,1);
      if (d > maxlen) {
        maxlen = d;
        kd = i;
      }
    }
    double cut = szs(kd,1) + maxlen/2.0;
    // Split
    for (int i = 0; i < n; i++) {
      double x = A(kd,is[i]);
      if (x <= cut) is1.push_back(is[i]);
      else is2.push_back(is[i]);
    }
  }
    
  // Calculate cluster data
  for (int i = 0; i < 2; i++) {
    ClusterData* cd = cdf->CreateClusterData();
    cd->Init(node->chldrn[i]->idxs, A, 0);
    node->chldrn[i]->cd = cd;
  }

  // Recurse
  for (int i = 0; i < 2; i++)
    Compute(node->chldrn[i]);
}

int ClusterTree::SetId (Node* node, int id) {
  node->id = id;
  for (int i = 0, n = node->chldrn.size(); i < n; i++)
    id = SetId(node->chldrn[i], id + 1);
  return id;
}

void ClusterTree::Leaves (vector<const Node*>& lvs) const {
  lvs.resize(0);
  Leaves(lvs, root);
}

void ClusterTree::Leaves (vector<const Node*>& lvs, Node* node) const {
  int n = node->chldrn.size();
  if (n == 0)
    lvs.push_back(node);
  else
    for (int i = 0; i < n; i++) Leaves(lvs, node->chldrn[i]);
}

void ClusterTree::Permutation (Vecui& p) const {
  int n = A.Size(2);
  p.resize(n);
  vector<const Node*> lvs;
  Leaves(lvs);
  for (UInt i = 0, k = 0; i < lvs.size(); i++)
    for (UInt j = 0; j < lvs[i]->idxs.size(); j++, k++)
      p[k] = lvs[i]->idxs[j];
}

// ---- HdCc

typedef vector<const Node*> ctnvec;

class HdCcCD : public ClusterData {
public:
  virtual void Init(const Vecui& idxs, const Matd& A, const Matd* wgt) = 0;

  unsigned int Nnodes () const { return nnodes; }

protected:
  unsigned int nnodes;
  Matd centroid;
  double radius;
};

class HdCcCDCentroid : public HdCcCD {
public:
  virtual void Init(const Vecui& idxs, const Matd& A, const Matd* wgt);

  void Centroid (Matd& c) const { c = centroid; }
  double Radius () const { return radius; }
};

class HdCcCDPointwise : public HdCcCDCentroid {
public:
  virtual void Init(const Vecui& idxs, const Matd& A, const Matd* wgt);

  // Diameter of this cluster
  double Diameter () const { return 2*radius; }
  // Shortest distance between a point in this cluster and one in cdp
  double Distance(const Hd* hd, const HdCcCDPointwise& cdp) const;

private:
  Vecui idxs;
  const Matd* A; // We can keep a pointer
};

class HdCcCDFactory : public ClusterDataFactory {
public:
  HdCcCDFactory (HdCc::DistanceMethod dm) : dm(dm) {}

  virtual ClusterData* CreateClusterData();

private:
  HdCc::DistanceMethod dm;
};

ClusterData* HdCcCDFactory::CreateClusterData () {
  switch (dm) {
  case HdCc::distCentroid:
    return new HdCcCDCentroid();
    break;
  case HdCc::distPointwise:
  case HdCc::distAxisAligned:
  default:
    // for now:
    return new HdCcCDPointwise();
    break;
  }
}

void HdCcCDCentroid::Init (const Vecui& is, const Matd& A, const Matd* wgt) {
  int n = is.size(), nd = A.Size(1);
  nnodes = n;

  // Centroid
  centroid.Resize(nd);
  centroid.Bzero();
  for (int i = 0; i < n; i++)
    for (int j = 1; j <= nd; j++)
      centroid(j) += A(j,is[i]);
  centroid *= 1.0 / n;

  // Radius
  radius = 0.0;
  for (int i = 0; i < n; i++) {
    double dist2 = 0.0;
    for (int j = 1; j <= nd; j++) {
      double d = A(j,is[i]) - centroid(j);
      dist2 += d*d;
    }
    if (dist2 > radius) radius = dist2;
  }
  radius = sqrt(radius);
}

void HdCcCDPointwise::Init (const Vecui& is, const Matd& A, const Matd* wgt) {
  HdCcCDCentroid::Init(is, A, wgt);
  idxs = is;
  this->A = &A;
}

// Minimum pointwise distance between two point clouds. Works only if the two
// point clouds do not overlap; otherwise, the correct value of 0 may not be
// returned.
double HdCcCDPointwise::
Distance (const Hd* hd, const HdCcCDPointwise& cdp) const {
  const Matd& A = *this->A;
  int n1 = idxs.size(), n2 = cdp.idxs.size(), nd = A.Size(1);

  //todo periodic
  double dist2 = -1.0;
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
#if 1
      double d2 = hd->Distance2(&A(1,idxs[i]), &A(1,cdp.idxs[j]), nd);
#else
      double d2 = 0.0;
      for (int k = 1; k <= nd; k++) {
        double a = A(k,idxs[i]) - A(k,cdp.idxs[j]);
        d2 += a*a; 
      }
#endif
      if (d2 < dist2 || dist2 < 0.0) dist2 = d2;
    }

  return sqrt(dist2);
}

Hd::iterator Hd::Begin () const { return blocks.begin(); }
Hd::iterator Hd::End () const { return blocks.end(); }
vector<hd::Block>::size_type Hd::NbrBlocks () const { return blocks.size(); }

void Hd::SetPeriodicBoundary (int dim, double lo, double hi) {
  if (dim > (int) _pbounds.Size(2)) {
    _pbounds.Reshape(2, dim);
    _is_periodic.resize(dim, false);
    _pwidths.resize(dim, 0.0);
  }
  _is_periodic[dim - 1] = true;
  _pbounds(1, dim) = lo;
  _pbounds(2, dim) = hi;
  _pwidths[dim - 1] = hi - lo;
}

static inline double amb2 (double a, double b) {
  double d = a - b; return d*d;
}

double Hd::Distance2 (const double* p1, const double* p2, int nd) const {
  double dist2 = 0.0;
  for (int i = 0; i < nd; i++) {
    double d2 = amb2(p1[i], p2[i]);
    if (i < (int) _is_periodic.size() && _is_periodic[i]) {
      double d2p = amb2(p1[i],
                        ((p1[i] < p2[i]) ?
                         p2[i] - _pwidths[i] :
                         p2[i] + _pwidths[i]));
      if (d2p < d2) d2 = d2p;
    }
    dist2 += d2;
  }
  return dist2;
}

HdCc::HdCc ()
  : eta(3.0), split_method(HdCc::splitMajorAxis),
    dist_method(HdCc::distCentroid)
{}

void HdCc::Opt_eta (double eta) throw (Exception) {
  if (eta <= 0) throw Exception("0 < eta");
  this->eta = eta;
}

bool HdCc::Admissible (const HdCcCD& n1cd, const HdCcCD& n2cd) const {
  if (n1cd.Nnodes() <= 1 || n2cd.Nnodes() <= 1) return false;

  switch (dist_method) {
  case distCentroid: {
    const HdCcCDCentroid& n1 = static_cast<const HdCcCDCentroid&>(n1cd);
    const HdCcCDCentroid& n2 = static_cast<const HdCcCDCentroid&>(n2cd);

    Matd c1(3), c2(3);
    n1.Centroid(c1);
    n2.Centroid(c2);
    int nd = c1.Size();

    double r1 = n1.Radius(), r2 = n2.Radius(),
      dist = Distance2(c1.GetPtr(), c2.GetPtr(), nd);
    dist = sqrt(dist) - (r1 + r2);

    //analysis
    // Magic number for testing weak admissibility. The construction of the
    // cluster tree in the symmetric case assures the clusters are
    // non-overlapping, which means any two trivially match the weak
    // admissibility condition (which is that they don't overlap) if their
    // centroids are not identical. (In fact, if they are identical, that means
    // n1cd and n2cd are the same object.)
    if (eta >= 1e20 && dynamic_cast<const HdCcSym*>(this))
      return Distance2(c1.GetPtr(), c2.GetPtr(), nd) > 0;
    
    return 2.0*min(r1, r2) < eta*dist;
  } break;

  case distPointwise:
  case distAxisAligned:
  default: {
    const HdCcCDPointwise& n1 = static_cast<const HdCcCDPointwise&>(n1cd);
    const HdCcCDPointwise& n2 = static_cast<const HdCcCDPointwise&>(n2cd);

    return min(n1.Diameter(), n2.Diameter()) <= eta*n1.Distance(this, n2);
  } break;
  }
}

const ctnvec* AssignCtnvec (const Node* n, ctnvec& spare) {
  if (!n->Chldrn().empty())
    return &n->Chldrn();
  spare.push_back(n);
  return &spare;
}

HdCcSym::HdCcSym (const Matd& A, const Matd* wgt)
  : A(A)
{}

void HdCcSym::Compute () { Compute(false); }

BNode* HdCcSym::Compute (bool want_bct) {
  BNode* bn = 0;

  // Cluster tree
  HdCcCDFactory cdf(dist_method);
  ClusterTree ct(&cdf, split_method, &A);
  ct.Compute();

  // Permutation p of the matrix B such that the blocks are all contiguous in
  // B(p,p).
  ct.Permutation(permp);
  // Inverse permutation.
  InvPermutation(permp, permpi);

  // Hierarchical decomposition.
  if (want_bct) {
    bn = Compute(&ct.Root(), &ct.Root());
  } else {
    ctnvec ns1, ns2;
    ns1.push_back(&ct.Root());
    ns2.push_back(&ct.Root());
    Compute(ns1, ns2);
  }

  return bn;
}

void HdCcSym::Compute (const ctnvec& ns1, const ctnvec& ns2) {
  for (int i = 0, m1 = ns1.size(); i < m1; i++) {
    const Node* n1 = ns1[i];
    bool hc1 = !n1->Chldrn().empty();
    for (int j = 0, m2 = ns2.size(); j < m2; j++) {
      const Node* n2 = ns2[j];
      if (n1->Id() > n2->Id()) continue;
      if (Admissible(static_cast<const HdCcCD&>(n1->CD()),
                     static_cast<const HdCcCD&>(n2->CD()))) {
        // n2 is far enough away from n1
        blocks.push_back(hd::Block(permpi, permpi, n1->Idxs(), n2->Idxs()));
        blocks.push_back(hd::Block(permpi, permpi, n2->Idxs(), n1->Idxs()));
      } else {
        // n2 is too near n1
        bool hc2 = !n2->Chldrn().empty();
        if (hc1 || hc2) {
          const ctnvec *a1, *a2;
          ctnvec spare;
          a1 = AssignCtnvec(n1, spare);
          a2 = AssignCtnvec(n2, spare);
          Compute(*a1, *a2);
        } else {
          blocks.push_back(hd::Block(permpi, permpi, n1->Idxs(), n2->Idxs()));
          if (n1->Id() != n2->Id())
            blocks.push_back(hd::Block(permpi, permpi, n2->Idxs(), n1->Idxs()));
        }
      }
    }
  }
}

BNode* HdCcSym::Compute (const Node* n1, const Node* n2) {
  hd::Block b(permpi, permpi, n1->Idxs(), n2->Idxs());
  BNode* nd = new BNode(b.r0, b.m, b.c0, b.n);

  bool hc1 = !n1->Chldrn().empty(), hc2 = !n2->Chldrn().empty();
  if (Admissible(static_cast<const HdCcCD&>(n1->CD()),
                 static_cast<const HdCcCD&>(n2->CD())) ||
#if 0
      // This condition makes the resulting leaves the same (I'm almost
      // positive. Testing gives identical results, but I haven't actually
      // proved that the two algorithms give the same results. The matter is one
      // of academic, but not practical, interest restricted just to the
      // decomposition.) as in the original version of HdCcSym::Compute above.
      !(hc1 || hc2)) {
#else
    // This condition makes each block cluster tree node always have 0 or 4
    // children. This is what we want if we intend to obtain an LU factorization
    // using the resulting block cluster tree.
    !hc1 || !hc2) {
#endif
    blocks.push_back(b);
    nd->block_idx = nblocks;
    nblocks++;
  } else {
    if (hc1 && hc2) {
      nd->chld[0] = Compute(n1->Chldrn()[0], n2->Chldrn()[0]);
      nd->chld[1] = Compute(n1->Chldrn()[0], n2->Chldrn()[1]);
      nd->chld[2] = Compute(n1->Chldrn()[1], n2->Chldrn()[0]);
      nd->chld[3] = Compute(n1->Chldrn()[1], n2->Chldrn()[1]);
    }
    else if (hc1 && !hc2) {
      nd->chld[0] = Compute(n1->Chldrn()[0], n2);
      nd->chld[1] = Compute(n1->Chldrn()[1], n2);
    } else { // (!hc1 && hc2)
      nd->chld[0] = Compute(n1, n2->Chldrn()[0]);
      nd->chld[1] = Compute(n1, n2->Chldrn()[1]);
    }
  }
  return nd;
}

HdCcNonsym::HdCcNonsym (const Matd& D, const Matd& R, const Matd* wgt)
  : D(D), R(R)
{}

void HdCcNonsym::Compute () {
  // Cluster tree
  HdCcCDFactory cdf(dist_method);
  ClusterTree ctd(&cdf, split_method, &D);
  ctd.Compute(false);
  ClusterTree ctr(&cdf, split_method, &R);
  ctr.Compute(false);

  // Permutations p,q of the matrix B such that the blocks are all contiguous
  // in B(p,q)
  ctd.Permutation(permq);
  InvPermutation(permq, permqi);
  ctr.Permutation(permp);
  InvPermutation(permp, permpi);

  // Hierarchical decomposition
  ctnvec ns1, ns2;
  ns1.push_back(&ctr.Root());
  ns2.push_back(&ctd.Root());
  Compute(ns1, ns2);
}

void HdCcNonsym::Compute (const ctnvec& ns1, const ctnvec& ns2) {
  for (int i = 0, m1 = ns1.size(); i < m1; i++) {
    const Node* n1 = ns1[i];
    bool hc1 = !n1->Chldrn().empty();
    for (int j = 0, m2 = ns2.size(); j < m2; j++) {
      const Node* n2 = ns2[j];
      if (Admissible(static_cast<const HdCcCD&>(n1->CD()),
                     static_cast<const HdCcCD&>(n2->CD()))) {
        // n2 is far enough away from n1
        blocks.push_back(hd::Block(permpi, permqi, n1->Idxs(), n2->Idxs()));
      } else {
        // n2 is too near n1
        bool hc2 = !n2->Chldrn().empty();
        if (hc1 || hc2) {
          const ctnvec *a1, *a2;
          ctnvec spare;
          a1 = AssignCtnvec(n1, spare);
          a2 = AssignCtnvec(n2, spare);
          Compute(*a1, *a2);
        } else {
          blocks.push_back(hd::Block(permpi, permqi, n1->Idxs(), n2->Idxs()));
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Hd.hpp public functions.

Hd::~Hd () {}

static bool IsSame (const Matd& D, const Matd& R) {
  // If R is empty, then we assume R == D.
  if (R.Size() == 0) return true;
  if (D.Size(1) != R.Size(1) || D.Size(2) != R.Size(2)) return false;
  const double *pD = D.GetPtr(), *pR = R.GetPtr();
  if (pD == pR) return true;
  for (UInt i = 0, n = D.Size(); i < n; i++)
    if (pD[i] != pR[i]) return false;
  return true;
}

Hd* NewHd (const Matrix<double>& D, const Matrix<double>& R,
           const Matrix<double>* pb, double eta) {
  assert(D.Size(1) == 3);
  assert(R.Size() == 0 || R.Size(1) == 3);

  HdCc::SplitMethod sm = HdCc::splitMajorAxis;

  hmmvp::HdCc* hd;
  if (IsSame(D, R)) hd = new HdCcSym(D);
  else hd = new HdCcNonsym(D, R);

  hd->Opt_eta(eta);
  hd->Opt_split_method(sm);
  if (sm == HdCc::splitAxisAligned)
    hd->Opt_distance_method(HdCc::distAxisAligned);
  else
    hd->Opt_distance_method(HdCc::distCentroid);

  if (pb)
    for (int i = 1; i <= 3; i++) {
      double lo = (*pb)(1, i), hi = (*pb)(2, i);
      if (lo < hi)
        hd->SetPeriodicBoundary(i, lo, hi);
    }

  hd->Compute();

  return hd;
}

Hd* NewHd (const Matrix<double>& D, const Matrix<double>* pb, double eta) {
  return NewHd(D, Matrix<double>(), pb, eta);
}

Hd* NewHdAxisAligned (const Matrix<double>& D, const Matrix<double>& R,
                      const Matrix<double>* pb, double eta) {
  assert(D.Size(1) == 3);
  assert(R.Size() == 0 || R.Size(1) == 3);

  HdCc::SplitMethod sm = HdCc::splitAxisAligned;

  hmmvp::HdCc* hd;
  if (IsSame(D, R)) hd = new HdCcSym(D);
  else hd = new HdCcNonsym(D, R);

  hd->Opt_eta(eta);
  hd->Opt_split_method(sm);
  if (sm == HdCc::splitAxisAligned)
    hd->Opt_distance_method(HdCc::distAxisAligned);
  else
    hd->Opt_distance_method(HdCc::distCentroid);

  if (pb)
    for (int i = 1; i <= 3; i++) {
      double lo = (*pb)(1, i), hi = (*pb)(2, i);
      if (lo < hi)
        hd->SetPeriodicBoundary(i, lo, hi);
    }

  hd->Compute();

  return hd;
}

Hd* NewHdAxisAligned (const Matrix<double>& D, const Matrix<double>* pb,
                      double eta) {
  return NewHdAxisAligned(D, Matrix<double>(), pb, eta);
}

void WriteHd (const Hd* hd, const std::string& hd_filename)
  throw (FileException)
{
  FILE* fid = fopen(hd_filename.c_str(), "wb");
  if (!fid) throw FileException();

  Vecui permp, permq;
  UInt n;
  hd->Permutations(permp, permq);
  n = permp.size();
  write(&n, 1, fid);
  write(&permp[0], n, fid);
  n = permq.size();
  write(&n, 1, fid);
  write(&permq[0], n, fid);

  n = hd->NbrBlocks();
  vector<UInt> data(4*n);
  UInt k = 0;
  for (Hd::iterator it = hd->Begin(), end = hd->End();
       it != end; ++it, k += 4) {
    const hd::Block& b = *it;
    data[k    ] = b.r0;
    data[k + 1] = b.m;
    data[k + 2] = b.c0;
    data[k + 3] = b.n;
  }
  write(&n, 1, fid);
  write(&data[0], data.size(), fid);

  fclose(fid);
}

Hd* NewHd (const std::string& hd_filename)
  throw (Exception, FileException)
{
  FILE* fid = fopen(hd_filename.c_str(), "rb");
  if (!fid) throw FileException();

  HdData* hd = new HdData();
  try {
    UInt n;
    read(&n, 1, fid);
    hd->Permp().resize(n);
    read(&hd->Permp()[0], n, fid);
    read(&n, 1, fid);
    hd->Permq().resize(n);
    read(&hd->Permq()[0], n, fid);

    read(&n, 1, fid);
    vector<UInt> data(4*n);
    read(&data[0], 4*n, fid);
    fclose(fid);
    for (UInt i = 0, k = 0; i < n; i++, k += 4)
      hd->AppendBlock(hd::Block(data[k], data[k+1], data[k+2], data[k+3]));
  } catch (const Exception& e) {
    delete hd;
    throw e;
  }

  return hd;
}

void DeleteHd(Hd* hd) { delete hd; }

}
