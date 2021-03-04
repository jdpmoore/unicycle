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

#ifndef INCLUDE_HMMVP_HELP_TEXT
#define INCLUDE_HMMVP_HELP_TEXT
/* Help text for programs. I write the text in a plain text file and then run
   this routine:

    function WrapPrintHelp (fn)
      fid = fopen(fn, 'r');
      state = 0;
      while (true)
        ln = fgetl(fid);
        if (~ischar(ln)) break; end
        if (numel(ln) >= 10 && strcmp(ln(1:7), '#define'))
          if (state == 1) fprintf('"\n\n'); end
          fprintf('%s \\\n"', ln);
          state = 1;
        else fprintf('%s\\n\\\n', ln); end
      end
      fprintf('"\n\n');
      fclose(fid);
    end
*/

#define _hmmvp_Header_text_ \
"hmmvp: Software to form and apply Hierarchical Matrices\n\
Version 1.2\n\
A.M. Bradley ambrad@cs.stanford.edu\n\
CDFM Group, Geophysics, Stanford University\n\
\n\
\n\
"

#define _hmmvpbuild_PrintHelp_text_ \
"hmmvpbuild\n\
  This message.\n\
\n\
hmmvpbuild <key-value file>\n\
  Build an H-matrix.\n\
    The key-value file describes the problem. The easiest way to create a\n\
  key-value file is to use matlab/kvf.m. You can create one in C++ by using the\n\
  interface described in util/include/KeyValueFile.hpp and linking against\n\
  libhmmvp.\n\
    The key-value file must contain at least this field:\n\
      command [string]: Right now, there is just one command available:\n\
        compress .\n\
    Type 'hmmvp help <command>' for help on each command.\n\
\n\
hmmvpbuild help <command>\n\
\n\
\n\
"

#define _hmmvpbuild_compress_PrintHelp_text_ \
"compress: Build an H-matrix. The key-value file must contain the following\n\
additional fields:\n\
\n\
  write_hmat_filename [string]: Filename prefix to which to write H-matrix\n\
    data. The extension .hm is appended if it is not already present.\n\
\n\
  eta [real scalar, [3]]: eta controls how point clusters are formed and\n\
    paired. Let c1 and c2 be two point clusters. For clusters c, c1, c2 and\n\
    points p, p1, p2 define:\n\
      centroid(c): Unweighted centroid of the points\n\
      dist(p1, p2): Euclidean distance between p1 and p2\n\
      radius(c): max_{point p in cluster c} dist(p, centroid(c)).\n\
    Then two clusters are considered to be sufficiently distant if\n\
      min(radius(c1), radius(c2)) <\n\
        eta/2 [dist(centroid(c1), centroid(c2)) - (radius(c1) + radius(c2))].\n\
\n\
  tol [real scalar, [1e-6]]: The error tolerance for the H-matrix\n\
    approximation. See next.\n\
\n\
  err_method [['mrem-fro'], 'brem-fro', or 'mrem-abs']: There are two parts to\n\
    each error method:\n\
      1. MREM vs BREM and\n\
      2. (Abs)olute or relative to ||G||_F, the (Fro)benius norm of the matrix.\n\
    If 'abs', then the error bound of the approximation Ga to G has the form\n\
      ||Ga - G||_F <= tol;\n\
    if 'fro', then\n\
      ||Ga - G||_F <= tol ||G||_F.\n\
    Cal the right-hand-side quantity 'err'. Then each has the form\n\
      ||Ga - G||_F <= err.\n\
    MREM and BREM differ in how they achieve this bound by assigning different\n\
    error tolerances to each block Ga_i. brem-fro implements\n\
      ||Ga_i - G_i||_F <= tol ||G_i||_F.\n\
    mrem-(fro/abs) implements\n\
      ||Ga_i - G_i||_F  <= sqrt(m_i n_i) / sqrt(M N) err,\n\
    where G is M x N and G_i is m_i x n_i.\n\
      MREM is the more efficient method of MREM and BREM in the following\n\
    sense. It achieves the greater compression for a given tol. But BREM may be\n\
    preferable when you want something closer to an element-wise, rather than\n\
    norm-wise, relative error. (This issue is somewhat complicated.) I recommend\n\
    you stick with the default unless you have a reason to switch.\n\
\n\
  Bfro [real scalar; optional]: The user's estimate of the Frobenius norm of the\n\
    matrix. If one is not provided but is needed, then an internal method\n\
    estimates the norm. The internal method will fail if the matrix does not\n\
    have a strong diagonal.\n\
\n\
  use_hmat_filename [string; optional]: Use an old H-matrix possibly to speed up\n\
    making this new one. hmmvp decides whether it can actually use the old file,\n\
    so there is no harm in specifying one. The extension .hm can be omitted.\n\
\n\
  allow_overwrite [[0] or 1]: Allow an H-matrix file to be overwritten by this\n\
    new one?\n\
\n\
  do_extra_tasks_only [[0] or 1]: Execute only what is in DoExtraTasks*().\n\
\n\
  write_Bfro_only [[0] or 1]: Write ||B||_F, either the computed estimate or the\n\
    value passed in field 'Bfro', to the file write_hmat_filename with extension\n\
    .Bfroest. The format is a single number written as %%1.15e\\n. Then exit. If\n\
    this option is not used, no .Bfroest file is written because, in that case,\n\
    it is more accurate to compute ||B||_F from the H-matrix directly.\n\
\n\
  nthreads [[1], >=1]: Number of OpenMP threads, if OpenMP is available and\n\
    enabled.\n\
\n\
  greens_fn [string]: The kernel to use. These are as follows:\n\
\n\
    inverse-r: A simple kernel for testing and demonstration.\n\
      X [3xN real array]:\n\
      order [real scalar]:\n\
      delta [real scalar]:\n\
\n\
    DEV Add additional kernel descriptions here.\n\
\n\
\n\
"

#endif
