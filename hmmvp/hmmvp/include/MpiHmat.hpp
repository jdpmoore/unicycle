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

/* Routines to parallelize the full-index version of Hmat::Mvp(). We assume that
   MPI_Init has already been called.
     Basic usage:

       MPI_Init(&argc, &argv);
       hmmvp::MpiHmat<float>* hm;
       try {
         hm = new Hm::MpiHmat<float>(hmat_fn, nrhs);
       } catch (const Exception& e) {
         cout << e.GetMsg() << endl;
         MPI_Finalize();
         exit(-1);
       }
       hm->Mvp(x, y, nrhs);
       
       delete hm;
       MPI_Finalize();

     Compile with mpic++. Run with a command of the sort (in Unix):
       mpirun -np 4 ./a.out
*/

#ifndef INCLUDE_HMMVP_MPI_HMAT
#define INCLUDE_HMMVP_MPI_HMAT

#include "util/include/Defs.hpp"
#include "util/include/Exception.hpp"

namespace hmmvp {
  class Hmat;
}

namespace hmmvp {
  using namespace std;
  using namespace util;

  typedef long long Blint;
  
  template<typename real>
  class MpiHmat {
  public:
    // For a pure-MPI program, nthreads is 1. But for a hybrid MPI-OpenMP one,
    // nthreads should be set to >1.
    MpiHmat(const string& filename, UInt ncol = 1, UInt nthreads = 1)
      throw (FileException, Exception);
    ~MpiHmat();

    // Number of rows.
    Blint GetM() const;
    // Number of columns.
    Blint GetN() const;
    // Number of nonzeros in the H-matrix approximation.
    Blint GetNnz() const;

    // y = B*x, where x has ncol columns. y is allocated by the caller.
    //   x and y are the full vectors.
    //   If perms are on, then x needs to be non-NULL only on the root.
    //   y must be non-Null and have the size of the full problem. On exit, its
    // values are valid only on the root.
    //   If the permutation routines below are not used, then it's safe to set x
    // to NULL on all ranks but 0.
    void Mvp(const real* x, real* y, UInt ncol) const throw (Exception);

  public:
    // Usually these can be ignored.
    void TurnOnPermute ()  { do_perms = true;  }
    void TurnOffPermute () { do_perms = false; }
    const Blint* GetQ() const;
    const Blint* GetP() const;

  private:
    hmmvp::Hmat* hm;
    mutable vector<real> workr;
    bool am_root;
    Blint nnz;
    bool do_perms;
  };

}

#include "MpiHmat_inl.hpp"

#endif
