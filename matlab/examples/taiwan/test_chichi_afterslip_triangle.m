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
%    [evl.t,evl.y]=evl.ode45((0:0.00274:2),y0,options);
%
% AUTHOR: 
% Sylvain Barbot (March 30, 2015), Earth Observatory of Singapore

clear all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

s2y=60*60*24*365;
y2s=1./s2y;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% triangle dislocations EGF
earthModelGimbutas12=greens.gimbutas12(30e3,1/4);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rcv=geometry.triangleReceiver('./faults/hsu15_receiver',earthModelGimbutas12);

% rate-and-state coefficients
rcv.a=0*rcv.l+1e-2;
rcv.b=rcv.a-1e-2;

% reference velocity (m/yr)
rcv.Vo=0*rcv.l+3e1;

% confining pressure (MPa)
rcv.sigma=0*rcv.l+1000;

% constrain rake?
rcv.isRakeConstraint=true;

% define a broad direction of afterslip
afterslipHeading=-90;
rcv.rake=atan2(rcv.dv(:,1)*sind(afterslipHeading)+rcv.dv(:,2)*cosd(afterslipHeading), ...
    rcv.sv(:,1)*sind(afterslipHeading)+rcv.sv(:,2)*cosd(afterslipHeading))*180/pi;
    

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

evt={geometry.coseismicTriangle('faults/hsu15',0,earthModelGimbutas12)};


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S A N I T Y   C H E C K                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if true
    figure(1000);clf;set(gcf,'name','slip distribution & frictional properties');
    subplot(1,2,1);cla;hold on
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    sc=2e3;
    quiver3(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3), ...
        sc*(rcv.sv(:,1).*cosd(rcv.rake)+rcv.dv(:,1).*sind(rcv.rake)), ...
        sc*(rcv.sv(:,2).*cosd(rcv.rake)+rcv.dv(:,2).*sind(rcv.rake)), ...
        sc*(rcv.sv(:,3).*cosd(rcv.rake)+rcv.dv(:,3).*sind(rcv.rake)), 0,'k');
    h=colorbar('South');
    title(h,'slip (m)');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    axis equal
    
    
    subplot(1,2,2);cla;hold on
    rcv.plotPatch((rcv.b-rcv.a).*rcv.sigma);
    
    h=colorbar('South');
    title(h,'friction properties');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    axis equal
    
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%             G R E E N S   F U N C T I O N S          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl')
    % build stress kernels for integral equation
    evl=ode.ratestrengthening(src,rcv,evt);
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.xc,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.xc,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.rcv=rcv;
    %evl.evt=evt;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%    C O S E I S M I C   S T R E S S   C H A N G E     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


if true
    tss=evl.evt{1}.Fss*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
        evl.evt{1}.Fds*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    tsd=evl.evt{1}.Fsd*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
        evl.evt{1}.Fdd*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    
    stressDrop=-sum(tss.*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
                    tsd.*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake)))/ ...
                sum(evl.evt{1}.src.slip);
	centroidDepth=-sum(evl.evt{1}.src.slip.*evl.evt{1}.src.xc(:,3))/ ...
                   sum(evl.evt{1}.src.slip);
           
    fprintf('stress drop: %f MPa\n',stressDrop);
    fprintf('centroid depth: %f km\n',centroidDepth/1e3);
    
    m=max(max(abs(tss)),max(abs(tsd)))/2;
    
    sc=2e3;
    figure(2000);clf;set(gcf,'name','coseismic stress change');
    subplot(1,2,1);cla;hold on
    
    evt{1}.plotPatch(tss);
    evt{1}.plotUnitVectors(sc);
    
    h=colorbar('South');
    xlabel(h,'shear stress in the strike direction (MPa)');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*m);
    axis equal
    
    
    subplot(1,2,2);cla;hold on
    evt{1}.plotPatch(tsd);
    %evt{1}.plotUnitVectors(sc);
    quiver3(evt{1}.xc(:,1),evt{1}.xc(:,2),evt{1}.xc(:,3), ...
        sc*(evt{1}.sv(:,1).*cosd(evt{1}.rake)+evt{1}.dv(:,1).*sind(evt{1}.rake)), ...
        sc*(evt{1}.sv(:,2).*cosd(evt{1}.rake)+evt{1}.dv(:,2).*sind(evt{1}.rake)), ...
        sc*(evt{1}.sv(:,3).*cosd(evt{1}.rake)+evt{1}.dv(:,3).*sind(evt{1}.rake)), 0,'k');
    
    quiver3(evt{1}.xc(:,1),evt{1}.xc(:,2),evt{1}.xc(:,3), ...
        sc/1e1*(evt{1}.sv(:,1).*tss+evt{1}.dv(:,1).*tsd), ...
        sc/1e1*(evt{1}.sv(:,2).*tss+evt{1}.dv(:,2).*tsd), ...
        sc/1e1*(evt{1}.sv(:,3).*tss+evt{1}.dv(:,3).*tsd), 0,'r');
    
    h=colorbar('South');
    title(h,'shear stress in the dip direction (MPa)');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*m);
    axis equal
    
    fprintf('coseismic stress change: review friction properties before simulation.\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% initial condition
y0=zeros(1,rcv.N*evl.dgf);

% integration options
options=ode.odeset( ...
    'Refine',1, ...
    'RelTol',1e-1, ...
    'InitialStep',1e-7, ...
    'verbose',true, ...
    'MaxStep',0.00137);
[evl.t,evl.y]=evl.ode45(0:0.00274:2,y0,options);

%% export to Paraview

if false
    evl.exportvtp(evl.t,evl.y,'./vtk');
end

%% surface data simulation

% plot GPS time series
if true
    % GPS data (network.dat contains a list of stations with name and coordinates)
    gps=unicycle.manifold.gpsReceiver('./gps/chichi_gps_sites.dat',evl,3);
    
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


%% map

% plot GPS in map view
if false
    scale=2e4;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    tcv.plotPatch()
    src.plotPatch()
    for k=1:length(gps.stationName)
        %text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.us((k-1)*gps.vecsize+1,timeIndex),scale*gps.us((k-1)*gps.vecsize+2,timeIndex));
    end
    xlabel('east (m)'), ylabel('north (m)'), zlabel('up (m)')
    axis equal tight
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end


%% fault slip time series plot

if true
    index={10,fix(mean(find(rcv.a-rcv.b>0)))};
    component=2;
    colors='rgbk';
    
    figure(6);clf;set(gcf,'name','time series of patch slip')
    subplot(3,1,1);cla;
    if 0<src.N
        plot(evl.t,evl.t*src.slip(1),[':k']);
    end
    subplot(3,1,2);cla;
    subplot(3,1,3);cla;
    for k=1:length(index)
        pos = (index{k}-1)*evl.dgf;
        subplot(3,1,1);hold on
        plot(evl.t,evl.y(pos+component,:),'-o','MarkerSize',2,'Color',colors(k));
        %plot(t,y(pos+2,:),'b-+');
        xlabel('time (yr)');ylabel('fault slip (m)');box on;
        
        subplot(3,1,2);hold on
        Dt=diff(evl.t);
        v=abs(diff(evl.y(pos+component,:))./Dt)*y2s;
        plot(evl.t(2:end),log10(v),'-o','MarkerSize',2,'Color',colors(k));
        plot(evl.t,log10(mean(Dt.*v)/mean(Dt)))
        xlabel('time (yr)');ylabel('log10 of slip velocity (m/s)');box on;
        
        subplot(3,1,3);hold on
        plot(evl.t,log10(abs(evl.y(pos+3,:))),'-o','MarkerSize',2,'Color',colors(k));
        xlabel('time (yr)');ylabel('log10 of shear stress (MPa)');box on;
        
    end
    subplot(3,1,1);
    if ~isempty(evt)
        for p=1:3
            subplot(3,1,p);
            for k=1:length(evt)
                plot(evt{k}.t0*[1 1],get(gca,'ylim'),'k-')
            end
        end
    end
end

%% static properties 3d plot

if true
    figure(7);clf;
    set(gcf,'name','source and receiver fault geometry')
    hold on
    timeIndex=ceil(length(evl.t)/2);
    rcv.plotPatch(sqrt(evl.y(1:evl.dgf:end,timeIndex).^2+evl.y(2:evl.dgf:end,timeIndex).^2));
    rcv.plotSlipVectors(evl.y(1:evl.dgf:end,timeIndex)',evl.y(2:evl.dgf:end,timeIndex)',1e3)
    %rcv.plotpatch(G*rcv.l./((rcv.b-rcv.a).*rcv.sigma)/1e3)
    %rcv.plotpatch(rcv.a-rcv.b)
    h=colorbar();
    ylabel(h,'slip (m)')
    
    %title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',t(k),t(end)));
    axis equal tight;%view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([38 48]), xlabel('x (km)'), ylabel('y (km)')
end



%% slip evolution movie

if true
    % bounds
    xlim=[-40 50]*1e3;
    ylim=[-30 65]*1e3;
    component=2; % dip slip
    
    figure(8);clf;
    sc=2e3;
    for k=2:15:length(evl.t)
        figure(8);cla;hold on
        %toplot=log10((evl.y(component:evl.dgf:end,k)-evl.y(component:evl.dgf:end,k-1))/(evl.t(k)-evl.t(k-1)));
        toplot=evl.y(component:evl.dgf:end,k);
        rcv.plotPatch(toplot);
        rcv.plotSlipVectors(evl.y(1:evl.dgf:end,k)',evl.y(2:evl.dgf:end,k)',sc)
        
        h=colorbar();
        ylabel(h,'slip (m)')
        
        title(sprintf('dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;%view([-30,20]);
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-16 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(max(abs(toplot)),1e-9)));
        box on, view([38 48]), xlabel('x (km)'), ylabel('y (km)')
        
        pause(0.0125)
    end
end
