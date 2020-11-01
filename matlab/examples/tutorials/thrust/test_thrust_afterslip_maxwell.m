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

s2y=60*60*24*365;
y2s=1./s2y;
G=30e3;
nu=1/4;

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);
minmax=@(x) [min(x),max(x)];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

flt=geometry.receiver('./faults/receiver.seg',greens.okada92(G,nu));

% rate-and-state coefficients
flt.a=0*flt.L+1e-2;
flt.b=flt.a-4e-3;

% reference velocity (m/yr)
flt.Vo=0*flt.L+1e+1;

% characteristic slip distance
flt.l=0*flt.L+0.3;

% confining pressure (MPa)
flt.sigma=0*flt.L+100;

% prevent afterslip in coseismic regions
flt.isRakeConstraint=true;
flt.Vrake=0*flt.L+90;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=[];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={geometry.coseismicPatch('faults/coseismic.flt',0,greens.okada92(G,nu))};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S H E A R   Z O N E S                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

shz=geometry.shearZoneReceiver('./faults/lowercrust',greens.shearZone16(G,nu));

shz.etaM=zeros(shz.N,1)+3.171e3; % MPa yr

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           G E O M E T R Y   C H E C K                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    subplot(1,2,1);cla;hold on
    flt.plotPatch(flt.a-flt.b);
    shz.plotShearZoneEdges();
    
    h=colorbar('South');
    ylabel(h,'Friction properties');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('(a-b) frictional properties')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(flt.b-flt.a)))
    axis equal
    
    subplot(1,2,2);cla;hold on
    shz.plotShearZoneEdges();
    flt.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    
    h=colorbar('South');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic slip distribution')
    set(gca,'view',[47 40]);
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl','var')
    tic
    % build stress kernels for integral equation
    evl=ode.rateStrengtheningBurgers(src,flt,shz,evt,'./kernel_afterslip_maxwell/');
    toc
else
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.flt=flt;
    evl.evt{1}.src=evt{1};
    evl.flt.dgf=5;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S T R E S S   C H E C K                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    subplot(1,2,1);cla;hold on
    % coseismic stress change
    toplot=evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    flt.plotPatch(toplot);
    shz.plotShearZoneEdges();
    
    h=colorbar('South');
    ylabel(h,'Friction properties');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic stress change in the dip direction')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/5e1)
    axis equal
    
    subplot(1,2,2);cla;hold on
    flt.plotPatch();
    
    toplot=evl.evt{1}.KL{2,6}*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    
    h=colorbar('South');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Stress change \sigma_{33}')
    set(gca,'view',[47 40]);
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                 S I M U L A T I O N                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% initial condition
y0=zeros(1,flt.N*evl.flt.dgf+shz.N*evl.shz.dgf);

tic
% integration options
options=ode.odeset('Refine',1,'RelTol',1e-13,'InitialStep',1e-5);
% creates variables evl.t and evl.y
% the following exports all computational time steps between t=0 and t=3;
%evl.ode45netcdf([0 3],y0,options);
% the following exports only at prescribed time steps 0,0.1,0.2,...,3
evl.ode45([0:0.1:3],y0,options);
toc
%evl.y=ncread('output.nc','yout');
%evl.t=ncread('output.nc','t');
%% export to Paraview

options.exportVtp=false;
if options.exportVtp
    evl.exportvtp(evl.t,evl.y,'./vtk');
end

%% surface data simulation

% plot GPS time series
if true
    % GPS data (network.dat contains a list of stations with name and coordinates)
    if ~exist('gps','var')
        gps=unicycle.manifold.gpsReceiver('./gps/dense_network.dat',evl,3);
    else
        gps.simulation(evl);
    end
    
    figure(1);clf;set(gcf,'Name','GPS time series')
    stationId=56;
    subplot(3,1,1);cla;
    hold on, box on;
    displacementComponent=1;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, east displacement',gps.stationName{stationId}));
    ylabel('East displacement (m)');
    
    subplot(3,1,2);cla;
    displacementComponent=2;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, north displacement',gps.stationName{stationId}));
    ylabel('North displacement (m)');
    
    subplot(3,1,3);cla;
    displacementComponent=3;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, vertical displacement',gps.stationName{stationId}));
    ylabel('Vertical displacement (m)'), xlabel('time (yr)');
end


%% GPS in map view

if true
    scale=5e5;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    subplot(2,1,1);cla;hold on; box on;
    timeIndex=ceil(length(evl.t));
    plot(gps.x(:,1),gps.x(:,2),'^');
    flt.plotPatch(sqrt(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2+evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2));
    flt.plotSlipVectors(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),1e5);
    shz.plotShearZoneEdges();
    h=colorbar();
    ylabel(h,'Slip (m)')
    
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
    title(sprintf('Displacements from afterslip at time %f yr',evl.t(timeIndex)));
    
    subplot(2,1,2);cla;hold on; box on;
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    flt.plotPatch();
    
    toplot=1e6*evl.y((evl.flt.N*evl.flt.dgf)+6:evl.shz.dgf:end,timeIndex);
    shz.plotShearZoneEdges(toplot);
    h=colorbar();
    ylabel(h,'Strain x10^{-6}')
    
    for k=1:length(gps.stationName)
        %text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.uv((k-1)*gps.vecsize+1,timeIndex),scale*gps.uv((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
    title(sprintf('Displacements from viscoelastic flow at time %f yr',evl.t(timeIndex)));
end

%% passive receiver faults for stress simulation

if true
    % passive receiver definition
    pflt=unicycle.geometry.passiveReceiver('./faults/passive_receiver',evl);
    
    figure(3);clf;set(gcf,'Name','passive receiver stress time series')
    subplot(2,1,1);cla;
    hold on, box on;
    plot(evl.t,pflt.tr(1,:));
    title('tension (MPa)');
    ylabel('tension (MPa)');
    
    subplot(2,1,2);cla;
    plot(evl.t,pflt.pr(1,:))
    title('pressure (MPa)');
    ylabel('pressure (MPa)');
    
end


%% time series plot

if true
    index={10,fix(mean(find(flt.a-flt.b>0)))};
    component=2;
    colors='rgbk';
    
    figure(6);clf;set(gcf,'name','Time series of patch slip')
    subplot(3,1,1);cla;
    if isobject(src)
        if 0<src.N
            plot(evl.t,evl.t*src.slip(1),[':k']);
        end
    end    
    subplot(3,1,2);cla;
    subplot(3,1,3);cla;
    for k=1:length(index)
        pos=(index{k}-1)*evl.flt.dgf;
        
        subplot(3,1,1);hold on
        plot(evl.t,evl.y(pos+component,:),'-o','MarkerSize',2,'Color',colors(k));
        
        xlabel('Time (yr)');ylabel('Fault slip (m)');box on;
        
        subplot(3,1,2);hold on
        Dt=diff(evl.t);
        v=abs(diff(evl.y(pos+component,:))./Dt)*y2s;
        plot(evl.t(2:end),log10(v),'-o','MarkerSize',2,'Color',colors(k));
        plot(evl.t,log10(mean(Dt.*v)/mean(Dt)))
        xlabel('Time (yr)');ylabel('log10 of slip velocity (m/s)');box on;
        
        subplot(3,1,3);hold on
        plot(evl.t,log10(abs(evl.y(pos+3,:))),'-o','MarkerSize',2,'Color',colors(k));
        xlabel('Time (yr)');ylabel('log10 of shear stress (MPa)');box on;
        
    end

    if ~isempty(evt)
        for p=1:3
            subplot(3,1,p);
            for k=1:length(evt)
                plot(evt{k}.t0*[1 1],get(gca,'ylim'),'k-')
            end
        end
    end
end

%% 3d fault slip plot

if true
    xlim=[-50 250]*1e3;
    ylim=[-200 180]*1e3;
    
    figure(7);clf;
    set(gcf,'name','Source and receiver fault geometry')
    timeIndex=ceil(length(evl.t)/2);
    
    subplot(2,1,1);gca;hold on;
    shz.plotShearZoneEdges();
    flt.plotPatch(sqrt(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2+evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2));
    flt.plotSlipVectors(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),4e4)
    
    h=colorbar();
    ylabel(h,'slip (m)')
    
    title(sprintf('Afterslip (m) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;%view([-30,20]);
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-80e3 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
    
    subplot(2,1,2);gca;hold on;
    
    toplot=evl.y((evl.flt.N*evl.flt.dgf)+15:evl.shz.dgf:end,timeIndex)*1e6;
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    %flt.plotPatch();
    
    h=colorbar();
    ylabel(h,'Strain (x10^{-6})')
    
    title(sprintf('Strain (x10^{-6}) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;%view([-30,20]);
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-80e3 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
end




%% slip evolution movie

if true
    % bounds
    xlim=[-50 250]*1e3;
    ylim=[-200 180]*1e3;
    component=2; % dip slip
    
    figure(8);clf;
    for k=2:15:length(evl.t)
        %figure(8);
        subplot(2,1,1);cla;hold on
        %toplot=log10((evl.y(component:evl.flt.dgf:end,k)-evl.y(component:evl.flt.dgf:end,k-1))/(evl.t(k)-evl.t(k-1)));
        toplot=evl.y(component:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k);
        flt.plotPatch(toplot);
        flt.plotSlipVectors(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k),evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k),1e4)
        shz.plotShearZoneEdges();
        h=colorbar();
        ylabel(h,'Slip (m)')
        
        title(sprintf('Dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-80 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
        
        subplot(2,1,2);cla;hold on
        %flt.plotPatch();
        toplot=evl.y(evl.flt.N*evl.flt.dgf+6:evl.shz.dgf:end,k)*1e6;
        shz.plotShearZoneEdges(toplot);
        shz.plotShearZoneEdges();
        h=colorbar();
        ylabel(h,'Strain (x10^{-6})')
        
        title(sprintf('Strain t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-80 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
        
        pause(0.0125)
    end
end
