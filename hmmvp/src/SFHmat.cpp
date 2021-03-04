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

#include <stdio.h>
#include <string>
#include "hmmvp/include/Hmat.hpp"
#include "hmmvp/include/SFHmat.h"

// Since C and Fortran can't handle C++ objects, we need to hide a global one in
// this file.
static hmmvp::Hmat* g_hm = 0;

int sfhmat_init_ (const char* hmat_fn, fint* ncol, fint* max_nthreads,
		  fint len_hmat_fn) {
  sfhmat_cleanup_();
  try {
    std::string my_hmat_fn(hmat_fn, len_hmat_fn);
    g_hm = hmmvp::NewHmat(my_hmat_fn, *ncol, *max_nthreads);
  } catch (const util::Exception& e) {
    fprintf(stderr, "sfhmat_init_: %s\n", e.GetMsg().c_str());
    g_hm = NULL;
    return -1;
  }
  return 0;
}

int sfhmat_get_size_ (fint* m, fint* n) {
  if (!g_hm) return -1;
  *m = g_hm->GetM();
  *n = g_hm->GetN();
  return 0;
}

int sfhmat_mvp_ (const Real* x, Real* y, fint* ncol) {
  if (!g_hm) return -1;
  try {
    g_hm->Mvp(x, y, *ncol);
  } catch (const util::Exception& e) {
    fprintf(stderr, "sfhmat_mvp_: %s\n", e.GetMsg().c_str());
    return -1;
  }
  return 0;
}

void sfhmat_cleanup_ () {
  if (g_hm) hmmvp::DeleteHmat(g_hm);
  g_hm = NULL;
}
