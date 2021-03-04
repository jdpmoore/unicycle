% DESCRIPTION:
% UNICYCLE example for afterslip simulation for the Sep 18, 2004 Mw 6
% Parkfield, CA earthquake.
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

options.sanityCheck=true;
if options.sanityCheck
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    rcv.plotPatch();
    %rcv.plotPatchIndex();
    for k=1:length(evt)
        fprintf('plotting slip distribution %s for time %f yr.\n',evt{k}.name,evt{k}.t0);
        evt{k}.plotPatch(evt{k}.slip);
    end
    h=colorbar(); box on
    ylabel(h,'coseismic slip (m)')
    
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('coseismic slip')
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

% initial condition
y0=zeros(1,rcv.N*evl.dgf);

% state variable
y0(1,5:evl.dgf:end)=0;

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5);
[evl.t,evl.y]=evl.ode45((0:0.0027:1),y0,options);

%% export to Paraview

options.exportVtp=false;
if options.exportVtp
    evl.exportvtp(evl.t,evl.y,'./vtk');
end

%% surface data simulation

options.plotGpsTimeSeries=true;
if options.plotGpsTimeSeries
    % GPS data (network.dat contains a list of stations with name and coordinates)
    gps=unicycle.manifold.okada85.gps('./gps/scign.dat',evl,3);
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
        scatter(gps.x(stationId,1),gps.x(stationId,2),2e2,gps.ur((stationId-1)*gps.vecsize+3,timeIndex)*1e3,'filled','MarkerEdgeColor','k');
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

%% profile and cross-section plots

options.plotHorizontalProfile=false;
if options.plotHorizontalProfile
    figure(3);clf;set(gcf,'Name','Horizontal profiles')
    hold on;box on
    depth=-mean(rcv.xc(:,3));
    evl.plotHorizontalProfiles(depth,1e-6,1e-4,0.25,2)
    %set(gca,'xlim',[60 120]*1e3)
    xlabel('along-strike distance (m)')
    ylabel('cumulative slip (m)')
end

options.plotVerticalProfile=false;
if options.plotVerticalProfile
    figure(4);clf;set(gcf,'Name','Vertical cross section')
    hold on;box on
    evl.plotVerticalProfiles(150e3,1e-6,1e-4,0.25,2)
    %set(gca,'ylim',[60 120]*1e3)
    xlabel('cumulative slip (m)')
    ylabel('down-dip distance (m)')
end

%% time series plot

options.plotTimeSeries=false;
if options.plotTimeSeries
    index={1,fix(mean(find(rcv.a>0)))};
    component=2;
    colors='rgbk';
    
    figure(5);clf;set(gcf,'name','time series of patch slip')
    subplot(4,1,1);cla;
    if 0<src.N
        plot(evl.t,evl.t*src.slip(1),[':k']);
    end
    subplot(4,1,2);cla;
    subplot(4,1,3);cla;
    subplot(4,1,4);cla;
    for k=1:length(index)
        pos = (index{k}-1)*evl.dgf;
        subplot(4,1,1);hold on
        plot(evl.t,evl.y(pos+component,:),'-o','MarkerSize',2,'Color',colors(k));
        %plot(t,y(pos+2,:),'b-+');
        xlabel('time (yr)');ylabel('fault slip (m)');box on;
        
        subplot(4,1,2);hold on
        Dt=diff(evl.t);
        v=abs(diff(evl.y(pos+component,:))./Dt)*y2s;
        plot(evl.t(2:end),log10(v),'-o','MarkerSize',2,'Color',colors(k));
        plot(evl.t,log10(mean(Dt.*v)/mean(Dt)))
        xlabel('time (yr)');ylabel('log10 of slip velocity (m/s)');box on;
        
        subplot(4,1,3);hold on
        plot(evl.t,log10(abs(evl.y(pos+3,:))),'-o','MarkerSize',2,'Color',colors(k));
        xlabel('time (yr)');ylabel('log10 of shear stress (MPa)');box on;
        
        subplot(4,1,4);hold on
        plot(evl.t,log10(evl.y(pos+5,:)),'-o','MarkerSize',2,'Color',colors(k));
        xlabel('time (yr)');ylabel('log10 of state variable (s)');box on;
    end
    subplot(4,1,1);
    if ~isempty(evt)
        for p=1:4
            subplot(4,1,p);
            for k=1:length(evt)
                plot(evt{k}.t0*[1 1],get(gca,'ylim'),'k-')
            end
        end
    end
end

%% static properties 3d plot

options.plotGeometry=true;
if options.plotGeometry
    figure(6);clf;
    set(gcf,'name','source and receiver fault geometry')
    hold on
    timeIndex=ceil(length(evl.t)/2);
    rcv.plotPatch(sqrt(evl.y(1:evl.dgf:end,timeIndex).^2+evl.y(2:evl.dgf:end,timeIndex).^2));
    %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,timeIndex)',evl.y(2:evl.dgf:end,timeIndex)',1e2)
    %rcv.plotpatch(G*rcv.l./((rcv.b-rcv.a).*rcv.sigma)/1e3)
    %rcv.plotpatch(rcv.a-rcv.b)
    h=colorbar();
    ylabel(h,'velocity (m/yr)')
    
    %title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',t(k),t(end)));
    axis equal tight; %view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
end



%% slip evolution movie

options.plotEvolution=true;
if options.plotEvolution
    % bounds
    xlim=[-30 10]*1e3;
    ylim=[-10 35]*1e3;
    sc=1e1;
    component=1; % strike slip
    
    figure(7);clf;
    for k=2:20:length(evl.t)
        figure(7);
        subplot(2,1,1);cla;hold on
        toplot=-evl.y(component:evl.dgf:end,k);
        rcv.plotPatch(toplot);
        %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,k)',evl.y(2:evl.dgf:end,k)',1e1);
        
        h=colorbar();
        ylabel(h,'afterslip (m)');
        title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-25 0]*1e3);
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([-40 28]), xlabel('x (m)'), ylabel('y (m)')
        
        subplot(2,1,2);cla
        toplot=log10((evl.y(component:evl.dgf:end,k)-evl.y(component:evl.dgf:end,k-1))/(evl.t(k)-evl.t(k-1))*y2s);
        rcv.plotPatch(toplot);
        
        h=colorbar();
        ylabel(h,'velocity in strike direction (m/s)');
        title(sprintf('velocity in dip direction (m/s) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-25 0]*1e3);
        set(gca,'clim',[-12 -1]);
        box on, view([-40 28]), xlabel('x (m)'), ylabel('y (m)')
        
        pause(0.0125)
    end
end
