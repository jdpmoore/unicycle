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
% Sylvain Barbot (March 30, 2015), Earth Observatory of Singapore

clear all
import unicycle.*

s2y=60*60*24*365;
y2s=1./s2y;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% homogeneous elastic half space, Poisson's solid
earthModelOkada92=greens.okada92(30e3,1/4);

% triangle dislocations EGF
earthModelGimbutas12=greens.nikkhoo15(30e3,1/4);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rcv=geometry.receiver('./faults/receiver.seg',earthModelOkada92);
tcv=rcv.toTriangleReceiver();
tcv.earthModel=earthModelGimbutas12;

% rate-and-state coefficients
tcv.a=0*tcv.l+1e-2;
tcv.b=tcv.a-4e-3;

index=find(1==omega((tcv.xc(:,2)-0e3)/100e3).*omega((tcv.xc(:,1)-88.6327e3)/75e3));
tcv.b(index)=tcv.a(index)+4e-3;

% reference velocity (m/yr)
tcv.Vo=0*tcv.l+1e+1;

% characteristic slip distance
tcv.l=0*tcv.l+0.3;

% loading direction
tcv.Vpl=0*tcv.l+0.03;

% loading direction
tcv.Vrake=0*tcv.l+90;

% confining pressure (MPa)
tcv.sigma=0*tcv.l+100;

tcv.observationPoints={ ...
    unicycle.geometry.observationPoint(226), ... % middle of velocity-weakening region
    unicycle.geometry.observationPoint(46) ...   % shallow velocity-strengthening patch
    }; 

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=geometry.source('./faults/thrust.flt',earthModelOkada92);

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
    tcv.plotPatch(tcv.b-tcv.a);
    tcv.plotPatchIndex();
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
    %evl=ode.rateandstatedamping(src,tcv,evt);
    evl=ode.rateandstate(src,tcv,evt,'./kernels_thrust_triangle/');
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    %assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    %assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.flt=tcv;
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
    cla;hold on
    toplot=evl.KK{2,2};
    toplot=toplot(144,:)+toplot(144+540,:);
    tcv.plotPatch(toplot);
    %tcv.plotPatchIndex();
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



%tcv.Vs=tcv.Vs*0+3e0*3.15e7;

% initial condition
y0=zeros(1,tcv.N*evl.flt.dgf);

% initial stress
%y0(1,4:evl.flt.dgf:end)=tcv.sigma.*(tcv.mu0+(tcv.a-tcv.b).*log(tcv.Vpl./tcv.Vo))+tcv.earthModel.G*tcv.Vpl./(2*tcv.Vs);

% state variable
%y0(1,5:evl.flt.dgf:end)=log(tcv.Vo./tcv.Vpl);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5);
% creates the variables evl.t and evl.y
evl.ode45([0 100],y0,options);

%% export to Paraview

if false
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
    ylabel('Vertical displacement (m)'), xlabel('Time (yr)');
end


%% map

% plot GPS in map view
if true
    scale=2e4;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
    tcv.plotPatch()
    src.plotPatch()
    %toplot=sum(gps.KO{2,1},2);
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex));
        %quiver(gps.x(k,1),gps.x(k,2),scale*toplot((k-1)*gps.vecsize+1),scale*toplot((k-1)*gps.vecsize+2));
    end
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
    title(sprintf('Surface displacements at time %f yr',evl.t(timeIndex)));
end

%% passive receiver faults for stress simulation

if false
    % passive receiver definition
    prcv=unicycle.geometry.passiveReceiver('./faults/passive_receiver',evl);
    
    figure(3);clf;set(gcf,'Name','passive receiver stress time series')
    subplot(2,1,1);cla;
    hold on, box on;
    plot(evl.t,prcv.tr(1,:));
    title('tension (MPa)');
    ylabel('tension (MPa)');
    
    subplot(2,1,2);cla;
    plot(evl.t,prcv.pr(1,:))
    title('pressure (MPa)');
    ylabel('pressure (MPa)');
    
end


%% fault slip time series plot

if true
    index={10,fix(mean(find(tcv.b-tcv.a>0)))};
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

if true
    figure(7);clf;
    set(gcf,'name','source and receiver fault geometry')
    hold on
    timeIndex=ceil(length(evl.t)/2);
    src.plotPatch();
    %tcv.plotPatch(tcv.b-tcv.a);
    %tcv.plotPatch(sqrt(evl.y(1:evl.flt.dgf:end,timeIndex).^2+evl.y(2:evl.flt.dgf:end,timeIndex).^2));
    %tcv.plotSlipVectors(evl.y(1:evl.flt.dgf:end,timeIndex)',evl.y(2:evl.flt.dgf:end,timeIndex)',1e2)
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
        tcv.plotPatch(toplot);
        src.plotPatch();
        
        h=colorbar();
        ylabel(h,'slip (m)')
        
        title(sprintf('dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;%view([-30,20]);
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-165 0]*1e3);
        clim=get(gca,'clim');
        set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        box on, view([-41 30]), xlabel('x (km)'), ylabel('y (km)')
        
        pause(0.0125)
    end
end
