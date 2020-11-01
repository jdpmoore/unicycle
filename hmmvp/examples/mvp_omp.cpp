// Example usage of hmmvp's matrix-vector product. Serial or OpenMP.

#include <iostream>
#include <math.h>

#include "hmmvp/include/Hmat.hpp"

using namespace std;

int main (int argc, char** argv) {
  if (argc != 2) {
    cerr << argv[0] << " [H-matrix filename]" << endl;
    return -1;
  }

  const char* fn = argv[1];
  const int nrhs = 2, max_nthreads = 4;

  hmmvp::Hmat* hm;
  try {
    hm = hmmvp::NewHmat(fn, nrhs, max_nthreads);
  } catch (const util::Exception& e) {
    cout << e.GetMsg() << endl;
    exit(-1);
  }

  const hmmvp::Blint n = hm->GetN(), m = hm->GetM();
  cout << "B is " << m << " by " << n << endl;
  cout << "||B||_F = " << sqrt(hm->NormFrobenius2()) << endl;

  // Make x vectors.
  double* x = new double[n*nrhs];
  for (size_t i = 0, N = n*nrhs; i < N; i++)
    x[i] = drand48() - 0.5;

  // Make y vector.
  double* y = new double[m*nrhs];

  // Compute y = A*x.
  hm->Mvp(x, y, nrhs);

  // I don't actually do anything with y here.

  delete[] y;
  delete[] x;
  hmmvp::DeleteHmat(hm);
}
