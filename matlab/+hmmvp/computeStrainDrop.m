function depsilon=computeStrainDrop(sourceFilename)
% function depsilon = COMPUTESTRAINDROP(sourceFilename) computes the strain drop 
% depsilon from a coseismic slip distribution using the expression
%
%         /           
%         | n . E . s dA
%         /           
%   dE  = --------------
%             /
%             | s dA
%             /
%
% where n is the fault normal vector, E is the strain tensor, s is the
% coseismic slip on the fault.
%
% example:
%
%    hmmvp.computeStrainDrop('faults/hill+14_0089.flt')
%
% AUTHORS:
% Priyamvada Nanjundiah (March 25, 2015)
% Sylvain Barbot (March 25, 2015)
%
% SEE ALSO: unicycle, hmmvp

import unicycle.geometry.source
import hmmvp.*

% earth model
earthModel=unicycle.greens.okada92(1/2,0);

% source and receiver faults
flt=source(sourceFilename,earthModel);

% hmmvp stress calculations
tol=1e-5;
nthreads=12;
base_fn='hmmvp/srcrcv/hm_';
assert(7==exist('hmmvp/srcrcv','dir'),'unicycle.hmmvp::error','temporary work space not found.');
h=hmmvp.unihmmvp('buildGreensFunctionsOkada92',flt,flt,1/2,0,base_fn,tol,nthreads);

% receiver stress interactions
[Kss.id,Kss.nnz]=hmmvp.hmmvp('init',h.Kss.hm_filename,nthreads);
[Ksd.id,Ksd.nnz]=hmmvp.hmmvp('init',h.Ksd.hm_filename,nthreads);
[Kds.id,Kds.nnz]=hmmvp.hmmvp('init',h.Kds.hm_filename,nthreads);
[Kdd.id,Kdd.nnz]=hmmvp.hmmvp('init',h.Kdd.hm_filename,nthreads);

rs=1:flt.N;
Kss.A=hmmvp.hmmvp('extract',Kss.id,rs,rs);
Ksd.A=hmmvp.hmmvp('extract',Ksd.id,rs,rs);
Kds.A=hmmvp.hmmvp('extract',Kds.id,rs,rs);
Kdd.A=hmmvp.hmmvp('extract',Kdd.id,rs,rs);

% stretch nhat . E . shat
nDotEStrike=Kss.A*(cosd(flt.rake).*flt.slip)+Kds.A*(sind(flt.rake).*flt.slip);
nDotEDip   =Ksd.A*(cosd(flt.rake).*flt.slip)+Kdd.A*(sind(flt.rake).*flt.slip);
stretch=nDotEStrike.*cosd(flt.rake)+nDotEDip.*sind(flt.rake);

% slip-averaged stress drop
depsilon=-sum(stretch.*flt.slip)/sum(flt.slip);

end
