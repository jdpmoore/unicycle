% program UNICYCLE solves the governing equations for fault
% slip evolution using rate-and-state friction with quasi-
% dynamics.
%
% friction equilibrium is enforced on the "receiver" faults,
% including the stress interaction between patches, stress 
% loading from "source" faults and sudden perturbations from 
% prescribed "coseismic" events.
%
% slip patches are defined in terms of position, orientation,
% and slip, as illustrated below:
%
%                 N (x1)
%                /
%               /| Strike
%   x1,x2,x3 ->@------------------------      (x2)
%              |\        p .            \ W
%              :-\      i .              \ i
%              |  \    l .                \ d
%              :90 \  S .                  \ t
%              |-Dip\  .                    \ h
%              :     \. | Rake               \
%              |      -------------------------
%              :             L e n g t h
%              Z (x3)
%
% input files for "source" and "receiver" faults are in the
% "./faults" directory. The convention (x1, x2, x3) for north,
% east and depth is for input files only. The coordinate system
% (x, y, z) for east, north and up is used for calculations and
% visualization internally.
%
% for afterslip studies, use empty source patches and add a
% coseismic perturbation:
%
%    src=geometry.source();
%    evt={geometry.coseismicpatch('./faults/coseismic.flt',0)};
%
% to obtain a afterslip distribution at given time steps (here daily 
% solutions for 2 years), use:
%
%    evl.ode45((0:0.00274:2),y0,options);
%
% AUTHOR: 
% Sylvain Barbot (May 24, 2013), Earth Observatory of Singapore

clear all
import unicycle.*

% elastic properties
G=30e3;
% Poisson's ratio (Poisson's solid)
nu=1/4;

% year sto seconds and seconds to years
s2y=60*60*24*365;
y2s=1./s2y;

% Radiation damping coefficient
rho=2.7e3;
dampingCoefficient=G/sqrt(G/rho)*y2s/50;



%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% homogeneous elastic half space, Poisson's solid
earthModel=greens.okada92(G,nu);



%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rcv=geometry.receiver('./faults/receiver.seg',earthModel);

% rate-and-state coefficients
rcv.a=0*rcv.L+1e-2;
rcv.b=rcv.a-4e-3;

index=find(1==omega((rcv.xc(:,2)-100e3)/100e3).*omega((rcv.xc(:,3)+20e3)/20e3));
rcv.b(index)=rcv.a(index)+4e-3;

% reference velocity (m/yr)
rcv.Vo=0*rcv.L+1e+1;

% characteristic slip distance
rcv.l=0*rcv.L+0.2;

% confining pressure (MPa)
rcv.sigma=0*rcv.L+240;

% confining pressure (MPa)
rcv.mu0=0*rcv.L+0.6;

fprintf('h* = %2.2e km\n',max(G*rcv.l./(rcv.b-rcv.a)./rcv.sigma/1e3))

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=geometry.source('./faults/transform.flt',earthModel);
src.slip=src.slip*10;
%src=geometry.source();

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={};%{geometry.coseismicpatch('./faults/coseismic.flt',0)};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

options.sanityCheck=false;
if options.sanityCheck
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    rcv.plotPatch(rcv.b-rcv.a);
    src.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch();
    end
    h=colorbar();
    ylabel(h,'friction properties');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('(b-a) frictional properties')
    set(gca,'view',[-54 52]);
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
    evl=ode.rateandstatedamping(src,rcv,evt,'./kernels_damping/');
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.rcv=rcv;
    evl.evt=evt;
    evl.previousVelocity=zeros(rcv.N,1);
    evl.dampingCoefficient=dampingCoefficient;
end

% initial condition
y0=zeros(1,rcv.N*evl.flt.dgf);

% initial shear stress mu*sigma
y0(1,3:evl.flt.dgf:end)=rcv.mu0.*rcv.sigma.*cosd(rcv.Vrake);
y0(1,4:evl.flt.dgf:end)=rcv.mu0.*rcv.sigma.*sind(rcv.Vrake);

% state variable (steady state)
y0(1,5:evl.flt.dgf:end)=0;

tic
% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5);
% creates variables evl.t and evl.y
evl.ode45([0 100],y0,options);
toc

%return

%% export to Paraview

options.exportVtp=false;
if options.exportVtp
    evl.exportvtp(evl.t,evl.y,'./vtk');
end

%% surface data simulation

options.plotGpsTimeSeries=true;
if options.plotGpsTimeSeries
    % GPS data (network.dat contains a list of stations with name and coordinates)
    gps=unicycle.manifold.gpsReceiver('./gps/dense_network.dat',evl,3);
    
    figure(1);clf;set(gcf,'Name','GPS time series')
    stationId=61;
    subplot(2,1,1);cla;
    hold on, box on;
    displacementComponent=1;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, east displacement',gps.stationName{stationId}));
    ylabel('east displacement (m)');
    
    subplot(2,1,2);cla;
    displacementComponent=2;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, north displacement',gps.stationName{stationId}));
    ylabel('north displacement (m)'), xlabel('time (yr)');
end


%% map

options.plotGpsMap=true;
if options.plotGpsMap
    scale=2e4;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    rcv.plotPatch()
    src.plotPatch()
    for k=1:length(gps.stationName)
        %text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('east (m)'), ylabel('north (m)'), zlabel('up (m)')
    axis equal tight
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end

%% profile and cross-section plots

options.plotHorizontalProfile=true;
if options.plotHorizontalProfile
    figure(3);clf;set(gcf,'Name','Horizontal profiles')
    hold on;box on
    depth=mean(rcv.xc(:,3));
    evl.plotHorizontalProfiles(depth,1e-2,3.17e-8,1,1)
    %set(gca,'xlim',[60 120]*1e3)
    xlabel('along-strike distance')
    ylabel('cumulative slip')
end

options.plotVerticalProfile=true;
if options.plotVerticalProfile
    figure(4);clf;set(gcf,'Name','Vertical cross section')
    hold on;box on
    dist=mean(rcv.xc(:,2));
    evl.plotVerticalProfiles(dist,1e-2,3.17e-8,1,1)
    %set(gca,'ylim',[60 120]*1e3)
    xlabel('cumulative slip')
    ylabel('down-dip distance')
end

%% time series plot

options.plotTimeSeries=true;
if options.plotTimeSeries
    index={10,fix(mean(find(rcv.b-rcv.a>0)))};
    component=1;
    colors='rgbk';
    
    figure(5);clf;set(gcf,'name','time series of patch slip')
    subplot(4,1,1);cla;
    if isobject(src)
        if 0<src.N
            plot(evl.t,evl.t*src.slip(1),[':k']);
        end
    end
    subplot(4,1,2);cla;
    subplot(4,1,3);cla;
    subplot(4,1,4);cla;
    for k=1:length(index)
        pos = (index{k}-1)*evl.flt.dgf;
        subplot(4,1,1);hold on
        plot(evl.t,evl.y(pos+component,:),'-o','MarkerSize',2,'Color',colors(k));
        %plot(t,y(:,pos+2),'b-+');
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
    src.plotPatch();
    rcv.plotPatch(sqrt(evl.y(1:evl.flt.dgf:end,timeIndex).^2+evl.y(2:evl.flt.dgf:end,timeIndex).^2));
    rcv.plotSlipVectors(evl.y(1:evl.flt.dgf:end,timeIndex),evl.y(2:evl.flt.dgf:end,timeIndex),1e2)
    %rcv.plotpatch(G*rcv.l./((rcv.b-rcv.a).*rcv.sigma)/1e3)
    %rcv.plotpatch(rcv.a-rcv.b)
    h=colorbar();
    ylabel(h,'slip (m)')
    
    %title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',t(k),t(end)));
    axis equal tight;%view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
end



%% slip evolution movie

options.plotEvolution=false;
if options.plotEvolution
    % bounds
    xlim=[-100 100]*1e3;
    ylim=[-50 250]*1e3;
    sc=1e1;
    
    figure(7);clf;
    for k=2:5:length(evl.t)
        figure(7);
        subplot(2,1,1);cla;hold on
        toplot=evl.y(1:evl.dgf:end,k);
        rcv.plotPatch(toplot);
        %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,k)',evl.y(2:evl.dgf:end,k)',1e2)
        %src.plotPatch();
        
        h=colorbar();
        ylabel(h,'slip (m)')
    
        title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;%view([-30,20]);
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-35 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[0 1]*(max(abs(toplot))+1e-9));
        box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
        
        subplot(2,1,2);cla
        toplot=log10((evl.y(1:evl.dgf:end,k)-evl.y(1:evl.dgf:end,k-1))/(evl.t(k)-evl.t(k-1)));
        rcv.plotPatch(toplot);
        %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,k)',evl.y(2:evl.dgf:end,k)',1e2)
        %src.plotPatch();
        
        h=colorbar();
        ylabel(h,'velocity (m/yr)')
    
        title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;%view([-30,20]);
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-35 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
        
        pause(0.0125)
    end
end
