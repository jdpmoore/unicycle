// Example usage of hmmvp's matrix-vector product. MPI.

#include <iostream>
#include <math.h>
#ifdef UTIL_OMP
#include <omp.h>
#endif

#include "hmmvp/include/Hmat.hpp"
#include "hmmvp/include/MpiHmat.hpp"

using namespace std;

// ||y1 - y2||_2 / ||y1||_2
template<typename T>
static double relerr (const T* y1, const T* y2, const size_t n) {
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < n; i++) {
    double d = y1[i] - y2[i];
    num += d*d;
    den += y1[i]*y1[i];
  }
  return sqrt(num/den);
}

int main (int argc, char** argv) {
  // Standard MPI startup code:
  MPI_Init(&argc, &argv);
  int nproc, pid;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  printf("nproc = %d  pid = %d\n", nproc, pid);

  if (argc != 2) {
    if (pid == 0) cerr << argv[0] << " [H-matrix filename]" << endl;
    MPI_Finalize();
    return -1;
  }

  const char* fn = argv[1];
  const int nrhs = 2;

  // Simple MPI H-matrix object.
  hmmvp::MpiHmat<double>* hm;
  try {
    // Test the hybrid MPI-OpenMP implementation if OpenMP is available.
    const int nthreads =
#ifdef UTIL_OMP
      omp_get_max_threads() > 1 ? 2 : 1
#else
      1
#endif
      ;
    hm = new hmmvp::MpiHmat<double>(fn, nrhs, nthreads);
  } catch (const util::Exception& e) {
    if (pid == 0) cerr << e.GetMsg() << endl;
    MPI_Finalize();
    return -1;
  }

  const hmmvp::Blint n = hm->GetN(), m = hm->GetM();
  if (pid == 0) cout << "B is " << m << " by " << n << endl;

  // Make x vector. Only the root has to have a valid one.
  double* x = NULL;
  if (pid == 0) {
    x = new double[n*nrhs];
    for (size_t i = 0, N = n*nrhs; i < N; i++)
      x[i] = drand48() - 0.5;
  }

  // Make y vector. All ranks must allocate the full y.
  double* y = new double[m*nrhs];

  // Compute y = A*x. Only the root has a valid y after this call.
  for (size_t ir = 0; ir < 1; ir++)
    hm->Mvp(x, y, nrhs);

  delete hm;

  // Compare result with the serial version.
  if (pid == 0) {
    do {
      // Make a serial (same as OpenMP) H-matrix object.
      hmmvp::Hmat* shm;
      try {
        shm = hmmvp::NewHmat(fn, nrhs);
      } catch (const util::Exception& e) {
        break;
      }

      double* sy = new double[m*nrhs];
      shm->Mvp(x, sy, nrhs);
      // Report the norm-wise relative difference between the two versions of y.
      cout << "serial vs MPI relerr: " << relerr(sy, y, n*nrhs) << endl;
      delete[] sy;
      hmmvp::DeleteHmat(shm);
    } while (0);
  }

  delete[] y;
  if (x) delete[] x;
  MPI_Finalize();
}
