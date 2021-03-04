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

#include "hmmvp/include/Hmat.hpp"

#include <iostream>
using namespace std;

typedef void* CHmat;

extern "C" {
int chmat_init (const char* hmat_fn, const int nrhs, const int max_nthreads,
                CHmat* chm) {
  hmmvp::Hmat* hm;
  try {
    hm = hmmvp::NewHmat(hmat_fn, nrhs, max_nthreads);
  } catch (const util::Exception& e) {
    cout << e.GetMsg() << endl;
    *chm = NULL;
    return -1;
  }
  *chm = static_cast<CHmat>(hm);
  return 0;
}

int chmat_get_info (const CHmat chm, int* m, int* n) {
  if (!chm) return -1;
  const hmmvp::Hmat* hm = static_cast<const hmmvp::Hmat*>(chm);
  *m = hm->GetM();
  *n = hm->GetN();
  return 0;
}

int chmat_mvp (CHmat chm, const double* x, double* y, const int ncol) {
  if (!chm) return -1;
  hmmvp::Hmat* hm = static_cast<hmmvp::Hmat*>(chm);
  hm->Mvp(x, y, ncol);
  return 0;
}

void chmat_cleanup (CHmat chm) {
  if (!chm) return;
  hmmvp::DeleteHmat(static_cast<hmmvp::Hmat*>(chm));
}
}
