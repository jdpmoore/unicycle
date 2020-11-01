function varargout = ex (varargin)
%DIR Example usage of hmmvp, driven by Matlab.
%   View this file and follow directions prefixed with %DIR in the order they
% are written.
%   Additionally, type 'help hmmvp' to get detailed information on all hmmvp
% functionality available in the Matlab interface.
%
%DIR In Matlab,
% >> addpath examples; addpath matlab;
  [varargout{1:nargout}] = feval(varargin{:});
end

%DIR To proceed, you need to build bin/hmmvpbuild.

%DIR On the Matlab command line, run
% >> b = ex('build');
function b = build ()
% Build a (very small) H-matrix.
  
  %DIR Modify the directory in which to place the H-matrix. The other options
  % can remain the same, though you might modify N, the (approximate) size of
  % the H-matrix.
  b.o = struct('dir','./tmp/',...
               'err_method', 'mrem-fro',...
               'tol', 1e-3,...
               'order', 3,...
               'geom', 'edgecube',...
               'N', 3000);

  % Make a toy problem.
  figure(1); clf;
  b.p = make_problem(b.o);
  
  % Write a key-value file describing the problem. b.c is a struct containing
  % the key-value pairs.
  b.c = write_kvf(b.p);
  
  %DIR In a shell, run one of these:
  write_shell_cmds(b.c.kvf);
end

%DIR To proceed, you need to build matlab/hmmvp. See exmvp_omp.cpp and
% exmvp_mpi.cpp to see how to compute MVP in C++ and exmvp_omp.c for C.

function analyze (b)
%DIR To analyze errors in a small H-matrix, in Matlab, run
% >> ex('analyze', b);
% EWRE = element-wise relative error
% NWRE = norm-wise relative error
  analyze_errors(b.p, b.c, 'MREM');
end

function demo_mvp (b)
%DIR To see how to compute a matrix-vector product (MVP) in Matlab, read the
% following code and run (with b from 'build')
% >> ex('demo_mvp', b);
  hm_filename = b.c.write_hmat_filename;
  % Load the H-matrix A into memory. id is a pointer to this H-matrix. Use 4
  % threads if they are available.
  id = hmmvp.hmmvp('init', hm_filename, 4);
  % Get the matrix's number of rows.
  m = hmmvp.hmmvp('getm', id);
  % Get the number of columns.
  n = hmmvp.hmmvp('getn', id);
  % Make a test vector.
  x = randn(n, 1);
  % Compute y = A*x.
  y = hmmvp.hmmvp('mvp', id, x);
  % y should be equal to size(A, 1).
  assert(numel(y) == m);
  % Clean up the memory for A. In practice, don't clean up until you're
  % completely done with this H-matrix.
  hmmvp.hmmvp('cleanup', id);
end

function demo_mvp_advanced (b)
%DIR Products of the form y(rs,:) = A(rs,cs)*x(cs,:) can be computed. Type
% >> ex('demo_mvp_advanced',b)
% to run this demo.
  hm_filename = b.c.write_hmat_filename;
  % Load the H-matrix A into memory. id is a pointer to this H-matrix. Use 4
  % threads if they are available. Additionally, allow up to 5 MVP in one
  % call.
  id = hmmvp.hmmvp('init', hm_filename, 4, 5);
  % Get the matrix's number of rows.
  m = hmmvp.hmmvp('getm', id);
  % Get the number of columns.
  n = hmmvp.hmmvp('getn', id);
 
  % Make a random set of vectors.
  X = randn(n, 5);
  % Do the basic full MVP:
  Y = hmmvp.hmmvp('mvp', id, X);

  % Now use just a subset of the rows:
  rs = 1:2:m;
  Y1 = hmmvp.hmmvp('mvp', id, X, rs);
  fprintf(['Should be the same to ~sqrt(eps) or eps, ', ...
           'depending on tolerance: %1.1e\n'],...
	  relerr(Y(rs,:), Y1(rs,:)));
  
  % Subset of cols:
  cs = 1:3:n;
  ncs = setdiff(1:n, cs);
  Y2a = hmmvp.hmmvp('mvp', id, X, [], cs);
  Y2b = hmmvp.hmmvp('mvp', id, X, [], ncs);
  fprintf('relerr: %1.1e\n', relerr(Y, Y2a + Y2b));
  
  % Subset of rows and cols.
  Y2a = hmmvp.hmmvp('mvp', id, X, rs, cs);
  Y2b = hmmvp.hmmvp('mvp', id, X, rs, ncs);
  fprintf('relerr: %1.1e\n', relerr(Y(rs,:), Y2a(rs,:) + Y2b(rs,:)));
  
  % For repeated calls to 'mvp' having the same row and column subsets, call
  %     hmmvp('ststate',id);
  % to save the internal state associated with the index sets. For small
  % index sets, this can offer quite a speedup. Call
  %     hmmvp('strelease',id);
  % when finished. NB: There is no error checking, so failing to release the
  % internal state will be a silent bug.
  hmmvp.hmmvp('ststate',id);
  for (i = 1:10)
    Y2a = hmmvp.hmmvp('mvp', id, X, rs, cs);
  end
  hmmvp.hmmvp('strelease',id);
  
  hmmvp.hmmvp('cleanup', id);
end

function b = recompress (b)
%DIR Recompress an H-matrix to increase its compression efficiency. Run
% >> ex('recompress', b);
  
  % Copy the old key-value struct to a new one, add a use_hmat_filename
  % field, and modify the output filenames.
  b.cr = b.c;
  b.cr.use_hmat_filename = b.cr.write_hmat_filename;
  b.cr.write_hmat_filename = [b.cr.write_hmat_filename(1:end-3) '_rc.hm'];
  b.cr.kvf = [b.cr.kvf(1:end-4) '_rc.kvf'];
  % Now make the tolerance looser. It can also be left the same. On larger
  % matrices, recompressing at the same tolerance can slightly improve
  % compression. For any tol >= the old one, recompression should be a lot
  % faster than the original compression.
  b.cr.tol = 1e2*b.cr.tol;
  % Write the new kvf-file.
  hmmvp.kvf('Write', b.cr.kvf, b.cr, 1);
  
  %DIR In a shell, run one of these:
  write_shell_cmds(b.cr.kvf);
end

% ------------------------------------------------------------------------------
% Private.

function p = make_problem (p)
% Make a test problem based on the simple kernel 1/sqrt(r^2 + delta)^order.
  if (nargin < 1) p = []; end
  p = popts(p, {{'dir' '.'};
                {'tol' 1e-5};    % Error tolerance on approx to B
                {'err_method' 'mrem-fro'};
                {'eta' 3};
                {'geom' 'edgecube'};
                {'N' 1000};      % B will be approximately NxN
                {'order' 3};     % Order of the singularity
                {'delta' 1e-4};  % For Plummer softening
                {'wiggle' 0};    % Wiggle the points a little?
                {'vary_charge' 0}}); % Vary the charges by this order of mag
  p.zero_diag = 0;
  if (p.delta == 0) p.zero_diag = 1; end
  
  switch (p.geom)
    case 'line'
      p.X = [linspace(-1,1,p.N); zeros(1,p.N); zeros(1,p.N)];
    case 'square'
      n = ceil(sqrt(p.N));
      p.N = n^2;
      x = linspace(-1,1,n);
      [X Y] = ndgrid(x,x);
      p.X = [X(:)'; Y(:)'; zeros(1,p.N)];
    case {'cube'}
      n = ceil(p.N^(1/3));
      p.N = n^3;
      x = linspace(-1,1,n);
      [X Y Z] = ndgrid(x,x,x);
      p.X = [X(:)'; Y(:)'; Z(:)'];
    case 'surfcube'
      n = ceil(sqrt(p.N/6));
      p.N = 6*n^2;
      dx = 2/n;
      x = -1+dx/2:dx:1-dx/2;
      o = ones(1,n^2);
      [X Y] = ndgrid(x,x);
      X = X(:)'; Y = Y(:)';
      p.X = [X -X X X o -o; Y Y o -o Y -Y; o -o Y -Y X X];
    case 'edgecube'
      n = ceil(p.N/12);
      p.N = 12*n;
      dx = 2/n;
      x = -1+dx/2:dx:1-dx/2;
      o = ones(1,n);
      p.X = [x x x x -o o -o o -o o -o o;
             -o o -o o x x x x -o -o o o;
             -o -o o o -o -o o o x x x x];
    otherwise
      error(sprintf('geom cannot be %s'),p.geom);
  end

  if (p.wiggle)
    p.X = p.X.*(1 + p.wiggle*randn(size(p.X)));
  end
  
  if (p.vary_charge)
    w = cos(2*pi*(1:p.N)/p.N);
    miw = min(w); maw = max(w);
    w = (w - miw)/(maw - miw);
    p.charge = 1 + w*(10^p.vary_charge - 1);
  else
    p.charge = ones(1, p.N);
  end

  plot3(p.X(1,:),p.X(2,:),p.X(3,:),'.');  
end

function c = write_kvf (p)
  bfn = sprintf(...
    '%sex_g%sN%dor%d_em%stol%deta%1.1f',...
    p.dir, p.geom(1:2), p.N, p.order, p.err_method(1:2), round(log10(p.tol)),...
    p.eta);
  c.command = 'compress';
  c.write_hmat_filename = [bfn '.hm'];
  c.kvf = [bfn '.kvf'];
  c.greens_fn = 'inverse-r';
  c = transfer_fields(c, p, {'err_method' 'tol' 'eta' 'X' 'order' 'delta'});
  c.allow_overwrite = 1;
  hmmvp.kvf('Write', c.kvf, c, 1);
end

function write_shell_cmds (kvf)
  variants = {'s' 'omp' 'mpi'};
  for (i = 1:2)
    fprintf('./bin/hmmvpbuild_%s %s\n', variants{i}, kvf);
  end
  fprintf('mpirun -n 4 ./bin/hmmvpbuild_mpi %s\n', kvf);
end

function analyze_errors (p, c, ttl)
  n = size(c.X, 2);
  if (n > 4000)
    fprintf('N is too large for analyze_errors.\n');
    return;
  end
  gf_new(p);
  G_true = gf_eval(1:n, 1:n);
  gf_cleanup();
  [id nnz] = hmmvp.hmmvp('init', c.write_hmat_filename);
  G = hmmvp.hmmvp('extract', id, 1:n, 1:n);
  hmmvp.hmmvp('cleanup', id);
  
  ewre = abs((G_true - G)./G_true);
  imagesc(log10(ewre)); colorbar;
  cf = n^2/nnz;
  title(sprintf(['%s: image is log_{10} EWRE\n',...
                 'NWRE: %1.2e (requested %1.1e)\n',...
                 'max EWRE: %1.2e\n',...
                 'compression %1.2fx'],...
                ttl, relerr(G_true, G), c.tol, max(ewre(:)), cf));
  set(gca, 'xtick', [], 'ytick', []);
end

function cf = calc_cf (hmfn)
  [id nnz] = hmmvp.hmmvp('init', hmfn);
  m = hmmvp.hmmvp('getm', id);
  n = hmmvp.hmmvp('getn', id);
  hmmvp.hmmvp('cleanup', id);
  cf = m*n/nnz;
end

function gf_new (p)
% Set up the Green's functions. p is from make_problem. The gf_* routines
% mimic the computations in GreensFnInverseR.cpp for testing in Matlab.
  global gf;
  gf.order = p.order;
  gf.delta = p.delta;
  gf.zero_diag = p.zero_diag;
  gf.X = p.X;
  gf.charge = p.charge;
end
  
function gf_cleanup ()
  clear global gf;
end
  
function B = gf_eval (rs, cs)
% Evaluate B(rs,cs) for the requested rows rs and columns cs.
  global gf;
  nr = length(rs);
  ns = length(cs);
  B = zeros(nr,ns);
  for (i = 1:ns)
    r2 = sum((repmat(gf.X(:,cs(i)),1,nr) - gf.X(:,rs)).^2);
    % For parameter eps in eps^2 = delta, this is called Plummer softening, at
    % least in gravity simulations.
    r = sqrt(r2 + gf.delta);
    if (gf.order > 0)
      B(:,i) = 1./r.^gf.order;
    else
      B(:,i) = log(r);
    end
    B(:,i) = B(:,i)*gf.charge(cs(i));
    if (gf.zero_diag)
      B(r2 == 0, i) = 0;
    end
  end
  if (~gf.zero_diag && gf.delta == 0)
    B(isinf(B)) = 0;
  end
end

function o = popt (o, fld, val)
  if (~isfield(o, fld)) o.(fld) = val; end
end

function o = popts (o, fvs)
  for (i = 1:numel(fvs))
    if (~isfield(o, fvs{i}{1})) o.(fvs{i}{1}) = fvs{i}{2}; end
  end
end

function sd = transfer_fields (sd, ss, flds)
  if (~iscell(flds)) flds = {flds}; end
  for (i = 1:numel(flds))
    if (isfield(ss, flds{i})) sd.(flds{i}) = ss.(flds{i}); end
  end
end

function re = relerr (a, b, p)
  if (nargin < 3) p = 'fro'; end
  re = norm(a - b, p)/norm(a, p);
  if (isnan(re))
    if (isempty(a) || all(a(:) == b(:))) re = 0; end
  end
end
