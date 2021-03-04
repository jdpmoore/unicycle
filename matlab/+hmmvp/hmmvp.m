
% hmmvp: Matrix-vector product and related operations with an H-matrix.
%
% [id nnz] = hmmvp('init',filename,[max_nthreads],[ncol])
%   Create the H-matrix B. Reat the data from filename.
%   If OpenMP is available, use at most max_nthreads when computing a
% matrix-vector product. max_nthreads defaults to 1.
%   ncol is the maximum number of columns in X when performing the operation
% Y = B*X.
%   id is used in subsequent calls to identify the H-matrix.
%   nnz is an optional output and is the number of numbers stored in the
% H-matrix. It can be compared with nnz(B), where B is the full BEM matrix.
%
% m = hmmvp('getm',id)
% n = hmmvp('getn',id)
%   Get the size of the matrix: [m n] = size(B).
%
% nthreads = hmmvp('threads',id,nthreads)
%   Use nthreads, if available. If nthreads > max_nthreads, it will be
% max_nthreads on exit. This function does not need to be called unless you
% want to make nthreads < max_nthreads.
%
% hmmvp('reorganizememory',id)
%   Reorganize the memory image of the H-matrix to accelerate the MVP. Once
% this function is called, subsequent calls to 'threads' triggers
% reorganization automatically.
%
% hmmvp('cleanup',id)
%   Clear the H-matrix from memory.
%
% clear hmmvp
%   Clear all H-matrices from memory.
%
% y = hmmvp('mvp',id,x,[rs],[cs])
%   Compute y = B*x.
%   rs,cs are optional subsets of rows, cols. indexing starts at 1. Just rs
% can be specified. Computes
%     y(sort(rs),:) = B(sort(rs),sort(cs))*x(sort(cs),:)           (1)
% or
%     y(sort(rs),:) = B(sort(rs),:)*x.                             (2)
% The full x vector must be provided (not x(cs,:)), and the full y vector
% (not y(sort(rs),:)) is returned. The entries of y not in rs are 0.
%
% y = hmmvp('tmvp',id,x,[rs],[cs])
%   Transpose MVP.
%
% hmmvp('stsave',id)
% hmmvp('strelease',id)
%   If a sequence of MVP of type (1) or (2) are computed using the same
% rs,cs, one can speed up the operation by saving certain state
% information. Use these functions as follows:
%     hmmvp('stsave',id);
%     % ... sequence of calls using the same rs,cs:
%     % ... y = hmmvp('mvp',id,x,rs,cs);
%     hmmvp('strelease',id);
%
% nf2 = hmmvp('fronorm2',id)
%   Compute the Frobenius norm (squared) of B.
%
% n1 = hmmvp('norm1',id)
%   Compute the 1-norm of B.
%
% Brc = hmmvp('extract',id,rs,cs).
%   Extract the block B(rs,cs) from the H-matrix. rs,cs need not be sorted,
% and results are ordered just as rs,cs are. Obviously rs,cs should be small,
% for this operation returns a matrix having length(rs)*length(cs) nonzero
% elements.
%
% hmmvp('pon',id)
% hmmvp('poff',id)
%   (Very specialized.) Optionally disable permutations associated with the
% H-matrix. Right now, this is appropriate only for AMB's iterative equation
% solver code.
%
% [I J S] = hmmvp('bblocks',id,[cutoff]);
%   Extract all the full-block blocks of the H-matrix and return a triple
% {I,J,S} that is used to create a sparse matrix. In Matlab, one creates
% the corresponding sparse matrix as follows:
%     A = sparse(I,J,S);
%   The optional parameter 'cutoff' is a value of the sort 1e-3 or 1e-2. Any
% values that are < cutoff*max(abs(S)) are dropped.
%
% SEE ALSO: hmmvp
