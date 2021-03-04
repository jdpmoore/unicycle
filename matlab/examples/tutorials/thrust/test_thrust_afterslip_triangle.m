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
minmax=@(x) [min(x),max(x)];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rflt=geometry.receiver('./faults/receiver.seg',greens.okada92(G,nu));
flt=rflt.toTriangleReceiver();
flt.earthModel=greens.nikkhoo15(G,nu);

% rate-and-state coefficients
flt.a=zeros(flt.N,1)+1e-2;
flt.b=flt.a-4e-3;

% reference velocity (m/yr)
flt.Vo=zeros(flt.N,1)+1e+1;

% characteristic slip distance
flt.l=zeros(flt.N,1)+0.3;

% confining pressure (MPa)
flt.sigma=zeros(flt.N,1)+100;

% prevent afterslip in coseismic regions
flt.isRakeConstraint=true;
flt.Vrake=zeros(flt.N,1)+90;

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
%           G E O M E T R Y   C H E C K                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','Geometry & frictional properties');
    subplot(1,2,1);cla;
    flt.plotPatch(flt.a-flt.b);
    
    h=colorbar('South');
    ylabel(h,'Friction properties');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Frictional properties (a-b)')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(flt.b-flt.a)))
    axis equal
    
    subplot(1,2,2);cla;
    flt.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    
    h=colorbar('South');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Slip')
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
    evl=ode.rateandstate(src,flt,evt,'./kernels_afterslip_triangle/');
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
    figure(1001);clf;set(gcf,'name','geometry & frictional properties');
    subplot(1,2,1);cla;
    % coseismic stress change
    toplot=sum(evl.KK{2,2},2);
    %toplot=toplot(284,:)+toplot(284+540,:);
    %toplot=evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    %toplot=evl.y(4:evl.flt.dgf:end,2);
    flt.plotPatch(toplot);
    
    h=colorbar('South');
    ylabel(h,'Stress (MPa)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic traction change in the dip direction')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))))
    axis equal
    
    subplot(1,2,2);cla;
    flt.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    
    h=colorbar('South');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic slip')
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
y0=zeros(1,flt.N*evl.flt.dgf);

% state variable
y0(1,5:evl.flt.dgf:end)=0;

tic
% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-9);
% creates variables evl.t and evl.y
% the following exports all computational time steps between t=0 and t=3;
evl.ode45([0 3],y0,options);
% the following exports only at prescribed time steps 0,0.1,0.2,...,3
%evl.ode45([0:0.1:3],y0,options);
toc

%% export to Paraview

options.exportVtp=false;
if options.exportVtp
    evl.exportvtp(evl.t,evl.y,'./vtk');
end

%% surface data simulation

% plot GPS time series
if true
    % GPS data (network.dat contains a list of stations with name and coordinates)
    gps=unicycle.manifold.gpsReceiver('./gps/dense_network.dat',evl,3);
    
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
    scale=4e5;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    flt.plotPatch(sqrt(evl.y(1:evl.flt.dgf:end,timeIndex).^2+evl.y(2:evl.flt.dgf:end,timeIndex).^2));
    flt.plotSlipVectors(evl.y(1:evl.flt.dgf:end,timeIndex),evl.y(2:evl.flt.dgf:end,timeIndex),2e4);
    h=colorbar();
    ylabel(h,'Slip (m)')
    
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('east (m)'), ylabel('north (m)'), zlabel('up (m)')
    axis equal tight
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end

%% passive receiver faults for stress simulation

if false
    % passive receiver definition
    prcv=unicycle.geometry.passiveReceiver('./faults/passive_receiver',evl);
    
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
    xlim=[-10 200]*1e3;
    ylim=[-200 200]*1e3;
    
    figure(7);clf;
    set(gcf,'name','Source and receiver fault geometry')
    hold on
    timeIndex=ceil(length(evl.t)/2);
    
    flt.plotPatch(sqrt(evl.y(1:evl.flt.dgf:end,timeIndex).^2+evl.y(2:evl.flt.dgf:end,timeIndex).^2));
    %flt.plotSlipVectors(evl.y(1:evl.flt.dgf:end,timeIndex)',evl.y(2:evl.flt.dgf:end,timeIndex)',4e4)
    
    h=colorbar();
    ylabel(h,'slip (m)')
    
    title(sprintf('Afterslip (m) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;%view([-30,20]);
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-35e3 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([53 56]), xlabel('x (km)'), ylabel('y (km)')
end



%% slip evolution movie

if true
    % bounds
    xlim=[-10 200]*1e3;
    ylim=[-200 200]*1e3;
    component=2; % dip slip
    
    figure(8);clf;
    for k=2:15:length(evl.t)
        figure(8);cla;hold on
        %toplot=log10((evl.y(component:evl.flt.dgf:end,k)-evl.y(component:evl.flt.dgf:end,k-1))/(evl.t(k)-evl.t(k-1)));
        toplot=evl.y(component:evl.flt.dgf:end,k);
        flt.plotPatch(toplot);
        flt.plotSlipVectors(evl.y(1:evl.flt.dgf:end,k),evl.y(2:evl.flt.dgf:end,k),4e4)
        
        h=colorbar();
        ylabel(h,'Slip (m)')
        
        title(sprintf('Dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-35 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([53 56]), xlabel('x (km)'), ylabel('y (km)')
        
        pause(0.0125)
    end
end
