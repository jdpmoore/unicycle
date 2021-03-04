function varargout = unihmmvp (varargin)
% Unicycle interface to hmmvp.
%
%   h = UNIHMMVP('buildGreensFunctionsOkada92',src,rcv,G,nu,base_fn,tol,nthreads)
%
% construct H-matrix approximations to the 8 matrices Kss, Ksd, Kds, Kdd,
% Fss, Fsd, Fds, Fdd.
%
% INPUT:
%   base_fn  - contains the directory and base filename for all output files.
%   tol      - controls the approximation accuracy. 1e-3 is very loose, 1e-5 is
%              about right, and 1e-7 gives very good accuracy. 
%              Compression is greater for bigger tol. 
%              Accuracy in the approximation A to B is such that
%              norm(A - B, 'fro') <= tol norm(B, 'fro').
%   nthreads - is the number of OpenMP threads to use when compressing the
%              matrices.
% 
%   unihmmvp('assess',h,evl)
%
% INPUT:
%   h        - is from build_hmatrices and evl is from ode.evolution. evl is optional. If
%              it is provided, assess norm-wise relative error of each matrix and the overall
%              error. The overall error should be no greater than tol. Whether or not evl is
%              provided, display the compression factors for each matrix and all together.
%
% add the following line to /Applications/MATLAB.app/bin/matlab
%
%   export PATH="$PATH:/Users/priyamva001/Documents/src/unicycle/hmmvp/bin"
%
% to run hmmvpbuild_omp.
%
% First version Feb 2014 AMB ambrad@cs.stanford.edu
%
% SEE ALSO: hmmvp, unicycle

  [varargout{1:nargout}] = feval(varargin{:});
  
end

function h = buildGreensFunctionsOkada92(src,rcv,G,nu,base_fn,tol,nthreads)
  if (nargin < 7), nthreads=1; end
  if (nargin < 6), tol=1e-6; end
  
  import hmmvp.*
  
  srcs='KF';
  comps='sd';
  first=1;
  Bfro=[];
  for (is=1:2)
    for (isc=1:2)
      for (irc=1:2)
        % The current matrix.
        fld = [srcs(is) comps(isc) comps(irc)];
        
        if (first)
          % We're using Kss to figure out the error control, so we need to
          % run it first.
          assert(strcmp(fld,'Kss'));
        end
        
        % Write the key-value files.
        kvf_fn=[base_fn '.kvf'];
        c.rcv_kvf_fn=[base_fn 'rcv.kvf'];
        c.src_kvf_fn=[base_fn 'src.kvf'];
        c.write_hmat_filename= ...
            [base_fn fld sprintf('tol%d', -round(log10(tol)))];
        c.tol=tol;
        c.nthreads=nthreads;
        if (~isempty(Bfro)) c.Bfro=Bfro; end
        c.allow_overwrite=1;
        c.greens_fn='okada92';
        c.mu=G;
        c.nu=nu;
        c.src_comp=comps(isc);
        c.rcv_comp=comps(irc);
        
        hmmvp.kvf('Write',c.rcv_kvf_fn,rcv,1);
        
        if (1==is)
          hmmvp.kvf('Write',c.src_kvf_fn,rcv,1);
        else
          hmmvp.kvf('Write',c.src_kvf_fn,src,1);
        end
        hmmvp.kvf('Write',kvf_fn,c,1);
        
        % Run hmmvpbuild.
        h.(fld).hm_filename=c.write_hmat_filename;
        if strcmp(srcs(is),'K') || 0<src.N
            [status,result]=system(sprintf('hmmvpbuild_omp %s;', kvf_fn));
            if 0~=status
                ME=MException('unicycle:runtimeError', ...
                    sprintf('%s\n\nerror running hmmvpbuild_omp %s',result,kvf_fn));
                throw(ME);
            end
        end
        
        if (first)
          % Get the Frobenius norm of Kss.
          id=hmmvp('init',h.(fld).hm_filename);
          Bfro=sqrt(hmmvp('fronorm2',id));
          hmmvp('cleanup',id);
          first=0;
        end
      end
    end
  end
end

function h = buildGreensFunctionsGimbutas12(src,rcv,G,nu,base_fn,tol,nthreads)
  if (nargin < 7), nthreads=1; end
  if (nargin < 6), tol=1e-6; end
  
  import hmmvp.*
  import unicycle.utils.*
  
  % options
  c.rfac=1;
  c.want_fullspace=0;
  c.tol=tol;
  c.nthreads=nthreads;
  c.allow_overwrite=1;
  c.greens_fn='gimbutas+12';
  c.mu=G; 
  c.nu=nu;
  % mesh points
  c.vertices=rcv.x';
  % triangle vertices
  c.triangles=rcv.vertices';
  
  labels={'s','d'};
  slipDirections={[1,0,0],[0,1,0]};
  tractionDirections={[1,0,0],[0,1,0]};
  
  textprogressbar('# hmmvp triangle stress kernels: ');
  for k=1:2
      for j=1:2
          textprogressbar(((k-1)*2+j-1)/4*100);
          
          fld=['K' labels{j} labels{k}];
          % dislocation in local coordinates
          % x is along strike, y is along dip, and z tensile
          c.disl=slipDirections{j};
          % receiver traction vector in local coordinates
          c.rcv=tractionDirections{k};
          % write key-value files
          c.write_hmat_filename=[base_fn fld sprintf('_tol%d', -round(log10(tol)))];
          c.write_hd_filename=[base_fn fld sprintf('_tol%d', -round(log10(tol)))];
          hmmvp.kvf('Write',[c.write_hmat_filename '.kvf'],c,1);
          % run hmmvpbuild
          h.(fld).hm_filename=c.write_hmat_filename;
          [status,result]=system(sprintf('hmmvpbuild_omp %s;', c.write_hmat_filename));
          if 0~=status
              ME=MException('unicycle:runtimeError', ...
                  sprintf('%s\n\nerror running hmmvpbuild_omp %s',result,c.write_hmat_filename));
              throw(ME);
          end
          % check that matrices can be accessed
          id=hmmvp.hmmvp('init',c.write_hmat_filename);
          hmmvp.hmmvp('cleanup', id);
      end
  end
  
  textprogressbar(100);
  textprogressbar('');
  
end



function assess (h, evl)
  if (nargin < 2) evl=[]; end

  err_tot=0;
  nrm_tot=0;
  h_bytes_tot=0;
  full_bytes_tot=0;
  flds=fieldnames(h);
  for (i=1:numel(flds))
    f=flds{i};

    d=dir(h.(f).hm_filename);
    h_bytes_tot=h_bytes_tot+d.bytes;

    id=hmmvp('init',h.(f).hm_filename);
    m=hmmvp('getm',id);
    n=hmmvp('getn',id);
    fnb=4*m*n; % 4 is for single-precision storage.
    cf=fnb/d.bytes;
    full_bytes_tot=full_bytes_tot+fnb;
    if (~isempty(evl))
      B=hmmvp('extract',id,1:size(evl.(f),1),1:size(evl.(f), 2));
    
      nrm_tot=nrm_tot+norm(B,'fro').^2;
      err_tot=err_tot+norm(B-evl.(f),'fro').^2;

      re=relerr(evl.(f), B);
      fprintf('%s %4dx%4d relerr %1.3e compression %1.2f\n',f,m,n,re,cf);
    else
      fprintf('%s %4dx%4d compression %1.2f\n',f,m,n,cf);
    end
    hmmvp('cleanup',id);
  end
  
  if (~isempty(evl))
    fprintf('full operator relerr %1.3e compression %1.2f\n', ...
            sqrt(err_tot/nrm_tot),full_bytes_tot/h_bytes_tot);
  else
    fprintf('full operator compression %1.2f\n',full_bytes_tot/h_bytes_tot);
  end
end

function re = relerr (a, b, p)
  if (nargin<3) p='fro'; end
  re=norm(a-b,p)/norm(a, p);
end

function fn = make_tmp_filename ()
% Make a temporary filename.
  gp=get_global_parms();
  k=0;
  %todo I could use random numbers so a linear search isn't needed.
  for (k=0:1000)
    fn=[gp.odir filesep 'tmp' num2str(k)];
    d=dir([fn '.*']);
    if (isempty(d)) return; end
  end
  error('Failed to find a tmp filename.');
end
