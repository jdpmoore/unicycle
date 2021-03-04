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

/* Very simple Fortran interface. */

#ifndef INCLUDE_HMMVP_SF_HMAT
#define INCLUDE_HMMVP_SF_HMAT

/* Adjust these to match the Fortran sizes. */
typedef double Real;
typedef long long int fint;

extern "C" {
  /* For gfortran, we use the following Fortran -> C conversion rules:
     1. All lowercase function name.
     2. Underscore at end of function name.
     3. If a char* is an argument, then there is an extra argument at the end
        that gives the string length.
     4. All arguments are pointers.

     Email me (ambrad@cs.stanford.edu) if this does not work with your Fortran
     compiler. Tell me the compiler you are using and I'll see what conventions
     must be used.  */
  int sfhmat_init_(const char* hmat_fn, fint* ncol, fint* max_nthreads,
		   /* Do not pass this next argument in your Fortran file; it
		      is automatically added by gfortran.  */
		   fint len_hmat_fn);
  int sfhmat_get_size_(fint* m, fint* n);
  int sfhmat_mvp_(const Real* x, Real* y, fint* ncol);
  void sfhmat_cleanup_();
}

#endif
