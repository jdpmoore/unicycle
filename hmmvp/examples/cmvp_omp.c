/* Example usage of hmmvp's matrix-vector product. Serial or OpenMP in C. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hmmvp/include/CHmat.h"

int main (int argc, char** argv) {
  const char *fn = argv[1];
  const int nrhs = 2, max_nthreads = 4;  
  CHmat chm;
  int m, n;
  size_t i, N;
  double *x, *y;

  if (argc != 2) {
    printf("%s [H-matrix filename]\n", argv[0]);
    return -1;
  }

  if (chmat_init(fn, nrhs, max_nthreads, &chm)) exit(-1);

  chmat_get_info(chm, &m, &n);
  printf("B is %d by %d\n", m, n);

  /* Make x vectors. */
  x = (double*) malloc(n*nrhs*sizeof(*x));
  for (i = 0, N = n*nrhs; i < N; i++)
    x[i] = drand48() - 0.5;

  /* Make y vector. */
  y = (double*) malloc(m*nrhs*sizeof(*y));

  /* Compute y = A*x. */
  chmat_mvp(chm, x, y, nrhs);

  /* I don't actually do anything with y here. */

  free(y);
  free(x);
  chmat_cleanup(chm);
}
