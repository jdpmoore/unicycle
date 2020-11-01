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

flt=geometry.receiver('./faults/receiver.seg',earthModel);

% rate-and-state coefficients
flt.a=0*flt.L+1e-2;
flt.b=flt.a-4e-3;

index=find(1==omega((flt.xc(:,2)-0e3)/100e3).*omega((flt.xc(:,1)-88.6327e3)/70e3));
flt.b(index)=flt.a(index)+4e-3;

% reference velocity (m/yr)
flt.Vo=0*flt.L+1e+1;

% characteristic slip distance
flt.l=0*flt.L+0.3;

% confining pressure (MPa)
flt.sigma=0*flt.L+100;

% loading rate (m/yr)
flt.Vpl=0*flt.L+0.03;

% loading rake (m/yr)
flt.Vrake=0*flt.L+90;

flt.observationPoints={ ...
    unicycle.geometry.observationPoint(226), ... % middle of velocity-weakening region
    unicycle.geometry.observationPoint(46) ...   % shallow velocity-strengthening patch
    }; 

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%src=geometry.source('./faults/thrust.flt',earthModel);
src=[];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    flt.plotPatch(flt.b-flt.a);
    src.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch();
    end
    h=colorbar();
    ylabel(h,'friction properties');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('(b-a) frictional properties')
    set(gca,'view',[47 40]);
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
if ~exist('evl','var')
    % build stress kernels for integral equation
    evl=ode.rateandstate(src,flt,evt,'./kernels_thrust/');
else
    % skip building stress kernel and update velocities & friction properties
    evl.flt=flt;
    evl.rcv=rcv;
    evl.evt=evt;
    evl.flt.dgf=5;
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         S T R E S S   I N T E R A C T I O N S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    toplot=evl.KK{1,1};
    toplot=toplot(144,:);
    flt.plotPatch(toplot);
    flt.plotPatchIndex();
    %src.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch();
    end
    h=colorbar();
    ylabel(h,'friction properties');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('(b-a) frictional properties')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/10);
    axis equal
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E V O L U T I O N                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%



% initial condition
y0=zeros(1,flt.N*evl.flt.dgf);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5);
% creates the variables evl.t and evl.y (export at computational time steps)
evl.ode45([0 600],y0,options);
% creates the variables evl.t and evl.y (export at prescribed time steps)
%evl.ode45([0:0.1:300],y0,options);

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
    ylabel('east displacement (m)');
    
    subplot(3,1,2);cla;
    displacementComponent=2;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, north displacement',gps.stationName{stationId}));
    ylabel('north displacement (m)');
    
    subplot(3,1,3);cla;
    displacementComponent=3;
    gps.plotTimeSeries(evl.t,stationId,displacementComponent);
    title(sprintf('station %s, vertical displacement',gps.stationName{stationId}));
    ylabel('vertical displacement (m)'), xlabel('time (yr)');
end


%% GPS in map view

if true
    scale=1e4;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    flt.plotPatch()
    if isobject(src)
        %src.plotPatch()
    end
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('east (m)'), ylabel('north (m)'), zlabel('up (m)')
    axis equal tight
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
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



%% profile and cross-section plots

options.plotHorizontalProfile=true;
if options.plotHorizontalProfile
    figure(4);clf;set(gcf,'Name','Horizontal profiles')
    hold on;box on ;
    depth=-mean(flt.xc(:,3));
    evl.plotHorizontalProfiles(depth,2e-8,1e-3,5,2)
    set(gca,'xlim',[60 120]*1e3)
    hold off; box off;
    xlabel('along-strike distance (m)')
    ylabel('cumulative slip (m)')
end

options.plotVerticalProfile=true;
if options.plotVerticalProfile
    figure(5);clf;set(gcf,'Name','Vertical cross section')
    hold on;box on
    evl.plotVerticalProfiles(148e3,2e-8,1e-3,5,2)
    set(gca,'ylim',[60 120]*1e3)
    xlabel('cumulative slip (m)')
    ylabel('down-dip distance (m)')
end

%% time series plot

if true
    index={10,fix(mean(find(flt.b-flt.a>0)))};
    component=2;
    colors='rgbk';
    
    figure(6);clf;set(gcf,'name','time series of patch slip')
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

if false
    figure(7);clf;
    set(gcf,'name','source and receiver fault geometry')
    hold on
    timeIndex=ceil(length(evl.t)/2);
    if isobject(src)
        src.plotPatch();
    end
    evl.flt.plotPatch(sqrt(evl.y(1:evl.flt.dgf:end,timeIndex).^2+evl.y(2:evl.flt.dgf:end,timeIndex).^2));
    evl.flt.plotSlipVectors(evl.y(1:evl.flt.dgf:end,timeIndex),evl.y(2:evl.flt.dgf:end,timeIndex),1e3)
    evl.flt.plotPatch(G*flt.l./((flt.b-flt.a).*flt.sigma)/1e3)
    evl.flt.plotPatch(flt.a-flt.b)
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

if true
    % bounds
    xlim=[-10 200]*1e3;
    ylim=[-200 200]*1e3;
    component=2; % dip slip
    
    figure(8);clf;
    for k=2:15:length(evl.t)
        figure(7);cla;hold on
        toplot=log10((evl.y(component:evl.flt.dgf:end,k)-evl.y(component:evl.flt.dgf:end,k-1))/(evl.t(k)-evl.t(k-1)));
        toplot=evl.y(component:evl.flt.dgf:end,k);
        flt.plotPatch(toplot);
        if isobject(src)
            %src.plotPatch();
        end
        h=colorbar();
        ylabel(h,'Slip (m)')
        
        title(sprintf('Dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;%view([-30,20]);
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-165 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
        
        pause(0.0125)
    end
end
