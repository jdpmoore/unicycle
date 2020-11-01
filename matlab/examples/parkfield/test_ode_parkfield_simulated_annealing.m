% DESCRIPTION:
% UNICYCLE example of using simulated annealing to optimize friction 
% parameters of faults for the case of the Sep 18, 2004 Mw 6 Parkfield, CA 
% earthquake.
%
% AUTHOR:
% Sylvain Barbot (September 15, 2013), Earth Observatory of Singapore

%clear all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

% homogeneous elastic properties
G=30e3;
nu=1/4; % Poisson's solid

s2y=60*60*24*365;
y2s=1./s2y;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rcv=geometry.receiver('./faults/receiver');

% rate-and-state coefficients
rcv.a=0*rcv.L+1e-2;
rcv.b=rcv.a-10e-3;

% reference velocity (m/yr)
rcv.Vo=0*rcv.L+1e1;

% characteristic slip distance
rcv.l=0*rcv.L+0.3;

% confining pressure (MPa)
rcv.sigma=0*rcv.L+400;

% is rake constrained?
rcv.isRakeConstraint=true;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=geometry.source();

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={geometry.coseismicpatch('./faults/barbot+12a.dat',0)};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

options.sanityCheck=false;
if options.sanityCheck
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    rcv.plotPatch();
    %rcv.plotPatchIndex();
    for k=1:length(evt)
        fprintf('plotting slip distribution %s for time %f yr.\n',evt{k}.name,evt{k}.t0);
        evt{k}.plotPatch(evt{k}.slip);
    end
    colorbar, box on
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('(b-a) frictional properties')
    set(gca,'view',[-36 24]);
    axis equal
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl')
    % build stress kernels for integral equation
    evl=ode.ratestrengthening(src,rcv,evt,G,nu);
    
    %
    % or try using the full rate-and-state equations
    %
    % evl=ode.rateandstate(src,rcv,evt,G,nu);
    %
    % and compare.
    
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.rcv=rcv;
    %evl.cap=2;
    evl.evt{1}.src.rake=evt{1}.rake;
    %evl.evt=evl.eventKernels(evt,rcv,G,nu);
end

if options.sanityCheck
    figure(1001);clf;set(gcf,'name','coseismic stress change');
    toplot=evl.evt{1}.Fss*evl.evt{1}.src.slip;
    rcv.plotPatch(toplot./((rcv.a-rcv.b).*rcv.sigma));
    %rcv.plotPatchIndex();
    h=colorbar(); box on;
    ylabel(h,'stress change (MPa)')
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[-36 24]);
    axis equal
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end

time=(0:0.0027:1);

gps=unicycle.manifold.okada85.gps('./gps/scign.dat',evl,2);
gps.d=gps.dataVectorFromTimeSeries('./gps/ps',time);

% initial condition
y0=zeros(1,rcv.N*evl.dgf);

% state variable
y0(1,5:evl.dgf:end)=0;

objectivefunction=@(m)(forwardmodel(m,evl,time,options,gps));

x0=[4e-3,1e0];
%objectivefunction(x0);
Mmax=2;
[x0,~]=unicycle.optim.sim_anl(objectivefunction,x0,[1e-3,0.001],[2e-2,3],Mmax)


%% surface data simulation

options.plotGpsTimeSeries=true;
if options.plotGpsTimeSeries
    % GPS data (network.dat contains a list of stations with name and coordinates)
    
    stationId=4;
    [gt,ge,gn,~,~]=textread(['./gps/ps/' gps.stationName{stationId}],...
        '%f %f %f %f %f','commentstyle','shell');
    
    figure(1);clf;set(gcf,'Name','GPS time series')
    
    subplot(3,1,1);cla;
    hold on, box on;
    displacementComponent=1;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    plot(gt-gt(1),ge,'k');
    set(gca,'xlim',[0 evl.t(end)]);
    plot(get(gca,'xlim'),[0 0],'k-');
    title(sprintf('station %s, east displacement',gps.stationName{stationId}));
    ylabel('east displacement (m)');
    
    subplot(3,1,2);cla;
    hold on, box on;
    displacementComponent=2;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    plot(gt-gt(1),gn,'k');
    set(gca,'xlim',[0 evl.t(end)]);
    plot(get(gca,'xlim'),[0 0],'k-');
    title(sprintf('station %s, north displacement',gps.stationName{stationId}));
    ylabel('north displacement (m)');
    
    subplot(3,1,3);cla;
    hold on, box on;
    displacementComponent=3;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    plot(get(gca,'xlim'),[0 0],'k-');
    title(sprintf('station %s, vertical displacement',gps.stationName{stationId}));
    ylabel('vertical displacement (m)'), xlabel('time (yr)');
end


%% map

options.plotGpsMap=true;
if options.plotGpsMap
    scale=2e5;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    time=evl.t(timeIndex);
    xlim=[-40 10]*1e3;
    ylim=[-10 30]*1e3;
    
    plot(gps.x(:,1),gps.x(:,2),'^');
    rcv.plotPatch()
    for stationId=1:length(gps.stationName)
        [gt,ge,gn,~,~]=textread(['./gps/ps/' gps.stationName{stationId}],...
        '%f %f %f %f %f','commentstyle','shell');
        [~,index]=min(abs(gt-gt(1)-time));
        text(gps.x(stationId,1)+2e3,gps.x(stationId,2)+2e3,gps.stationName{stationId});
        %scatter(gps.x(stationId,1),gps.x(stationId,2),2e2,gps.ur((stationId-1)*gps.vecsize+3,timeIndex)*1e3,'filled','MarkerEdgeColor','k');
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*gps.ur((stationId-1)*gps.vecsize+1,timeIndex),scale*gps.ur((stationId-1)*gps.vecsize+2,timeIndex));
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*ge(index),scale*gn(index),'k');
        plot.ellipse(scale*0.005,scale*0.005,0,gps.x(stationId,1)+scale*ge(index),gps.x(stationId,2)+scale*gn(index),'k');
    end
    xlabel('east (m)'), ylabel('north (m)')
    set(gca,'xlim',xlim,'ylim',ylim,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    h=colorbar();
    ylabel(h,'vertical displacement (mm)')
    %axis square
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end



