function dtau=computeStressDrop(sourceFilename,G,nu)
% function dtau = HMMVP.COMPUTESTRESSDROP(sourceFilename,G,nu) computes the stress drop
% dtau from a coseismic slip distribution using
%
%         /
%         | tau * s dA
%         /
%   dtau = --------------
%             /
%             | s dA
%             /
%
% where tau is the amplitude of shear stress on the fault, s is the
% coseismic slip on the fault.
%
% example:
%
%    hmmvp.computeStressDrop('faults/hill+14_0089.flt',30e3,0.25)
%
% AUTHOR:
% Sylvain Barbot (May 8, 2014), Earth Observatory of Singapore
%
% SEE ALSO: unicycle, hmmvp

import unicycle.geometry.source
import hmmvp.*

% earth model
earthModel=unicycle.greens.okada92(G,nu);

% source and receiver faults
flt=source(sourceFilename,earthModel);

% hmmvp stress calculations
tol=1e-5;
nthreads=12;
base_fn='hmmvp/srcrcv/hm_';
assert(7==exist('hmmvp/srcrcv','dir'),'unicycle.hmmvp::error','temporary work space not found.');
h=hmmvp.unihmmvp('buildGreensFunctionsOkada92',flt,flt,G,nu,base_fn,tol,nthreads);

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

% traction
ts=Kss.A*(cosd(flt.rake).*flt.slip)+Kds.A*(sind(flt.rake).*flt.slip);
td=Ksd.A*(cosd(flt.rake).*flt.slip)+Kdd.A*(sind(flt.rake).*flt.slip);
tau=(ts.*cosd(flt.rake)+td.*sind(flt.rake));

% slip-averaged stress drop
dtau=-sum(abs(flt.slip).*tau)/sum(abs(flt.slip));

end
