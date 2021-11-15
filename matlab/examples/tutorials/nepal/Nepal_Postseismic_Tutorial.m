% Script Nepal_Postseismic_Tutorial simulates postseismic deformation following the 
% 2015 Mw 7.4 Gorkha earthquake with afterslip on the Himalayan Main 
% Frontal Thrust using the 3D fault geometry of Hubbard et
% al. (2016) and viscoelastic flow in the lower crust. 
%
% Friction dynamic equilibrium is enforced on the "receiver" fault, which
% is broken down into triangular elements. The geometry is described in two 
% files,
%
%   ./faults/qiu+15_1_receiver.ned
%   ./faults/qiu+15_1_receiver.tri
%
% The qiu+15_1_receiver.ned file contains a list of points that can be
% used to define the vertices of the triangular elements. The 
% qiu+15_1_receiver.tri file contains the list of triangular elements
% each defined by the index of three vertices, as illustrated below:
%
%        P1 (x1,x2,x3)  @---------@   P2 (x1,x2,x3)
%                        \       /
%                         \     /
%                          \   /
%                           \ /
%                            @ 
%       
%                      P3 (x1,x2,x3)
% 
% The convention (x1, x2, x3) for north, east and depth is for input files. 
% The coordinate system (x, y, z) for east, north and up is used for
% internal calculations and visualization in Matlab.
%
% The equation of state for a power law Burger's rheology
% used to predict the viscoelastic relaxation of the lower crust. The
% geometry is described in the file
%
%   ./faults/lc+sz_lvl.shz
%
% that contains a list of shear zones described as follows:
%
%
%                      N (x1)
%                     /
%                    /| strike (theta)          E (x2)
%        q1,q2,q3 ->@--------------------------+
%                   |                        w |     +
%                   |                        i |    /
%                   |                        d |   / s
%                   |                        t |  / s
%                   |                        h | / e
%                   |                          |/ n
%                   +--------------------------+  k
%                   :       l e n g t h       /  c
%                   |                        /  i
%                   :                       /  h
%                   |                      /  t
%                   :                     /
%                   |                    +
%                   Z (x3)
%
%
% The loading rate is imposed on each patch so that 1) the norm of the
% local loading rate is uniform and 2) the horizontal direction of loading
% is uniform.
%
%   flt.Vpl=0.02*ones(flt.N,1);
%   azimuth=-160;
%   azimuthVector=[0*flt.l+sind(azimuth),0*flt.l+cosd(azimuth),0*flt.l];
%   flt.Vrake=atan2(sum(flt.dv.*azimuthVector,2),sum(flt.sv.*azimuthVector,2))*180/pi;
%
% The calculation is carried out over several years and the time
% steps for daily output are specified explicitely:
%
%    [evl.t,evl.y]=evl.ode45(0:0.002739:3,y0,options);
%
% AUTHORS: 
% James Moore (6th August 2017), Earth Observatory of Singapore

clear all
close all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

load './gmt/boundaries.dat';

gpsdata=cell(4,1);
gpsdata{1} = importdata('./gps/KKN4out.dat');
gpsdata{2} = importdata('./gps/NASTout.dat');
gpsdata{3} = importdata('./gps/CHLMout.dat');
gpsdata{4} = importdata('./gps/SNDLout.dat');

s2y=60*60*24*365;
y2s=1./s2y;
colors='kmcrgby';
G=30e3;
nu=1/4;
minmax=@(x) [min(x),max(x)];
heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

flt=geometry.triangleReceiver('./faults/qiu+15_1_receiver',greens.nikkhoo15(G,nu));

% rate-and-state coefficients
flt.a=0*flt.l+1e-2;
flt.b=flt.a-8e-3;

strike=-72;
s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-0e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)+0e3)*sind(strike);

% reference velocity (m/yr)ls
flt.Vo=0*flt.l+1e-6;

% characteristic slip distance
flt.l=0*flt.l+7.5e-2;

% confining pressure (MPa)
flt.sigma=0*flt.l+1500;

% plate loading
flt.Vpl=0*flt.l+0.02;
azimuth=-160;
azimuthVector=[0*flt.l+sind(azimuth),0*flt.l+cosd(azimuth),0*flt.l];
flt.Vrake=atan2(sum(flt.dv.*azimuthVector,2),sum(flt.sv.*azimuthVector,2))*180/pi;
flt.rake=flt.Vrake;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={geometry.coseismicTriangle('./faults/qiu+15_1',0,greens.nikkhoo15(G,nu))};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%      R E C E I V E R   S H E A R   Z O N E S         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

shz=geometry.shearZoneReceiver('./faults/lc+sz',greens.shearZone16(G,nu));

pos=1:shz.levels{1}.nShearZone;

s=+(shz.x(pos,2)+250e3)*cosd(strike)+(shz.x(pos,1)-0e3)*sind(strike);
d=+(shz.x(pos,1)+140e3)*cosd(strike)-(shz.x(pos,2)+0e3)*sind(strike);

ramp=@(x) x.*omega(x-1/2)+heaviside(x-1);
shz.x(pos,3)=shz.x(pos,3)-35e3*ramp(d/1e5);
shz.xc(pos,3)=shz.xc(pos,3)-35e3*ramp(d/1e5);

% rheology
shz.n=0*shz.n+1;
shz.m=0*shz.n+1;
shz.etaM=shz.etaM*0+1e20/(365*24*60*60)/1e6;
%shz.T0 = importdata('temps2.dat');
%shz.etaM=arr(shz.T0);
shz.etaK=1e-2*shz.etaM;

% plate loading
shz.epsilonPlate=[0,1,0,0,0,0]*1e-15;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            G E O M E T R Y   C H E C K               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if true
    figure(1000);clf;set(gcf,'name','geometry & properties');
    hold on
    
    flt.plotPatch();
    evt{1}.plotPatch(evt{1}.slip);
    
    shz.plotShearZoneEdges();
    shz.plotShearZoneEdges(log10(3.1536e+13*shz.etaK));
        
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    colormap jet
    h=colorbar('southoutside');
    ylabel(h,'Coseismic slip (m) / Viscosity (Pa s)');
    box on, grid on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('Afterslip and Viscosity')
        
    set(gca,'view',[-227 42]);
    %set(gca,'view',[69.6000   -0.4000]);
    axis equal
    set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3)
    set(gca,'clim',[0 1]*max(abs(get(gca,'clim'))));
   %set(gca,'clim',[14.5 24]);
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S T R E S S   K E R N E L S             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl','var')
    % build stress kernels for integral equation
    tic
    evl=unicycle.ode.rateStrengtheningPower([],flt,shz,evt,'./kernels/');
    toc
else
    % skip building stress kernel and update velocities & friction properties
    evl.flt=flt;
    evl.flt.dgf=4;
    evl.shz=shz;
    evl.shz.dgf=18;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%          S T R E S S   I N T E R A C T I O N         %
%               V I S U A L I Z A T I O N              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if true
    figure(1001);clf;set(gcf,'name','Stress interactions');
    subplot(2,1,1);gca;hold on
    slip=evt{1}.slip;
     toplot=evl.evt{1}.KK{1,2}*(slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(slip.*sind(evt{1}.rake));
   % toplot=evl.evt{1}.KK{1,3}*(slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,3}*(slip.*sind(evt{1}.rake));
    flt.plotPatch(toplot);

    plot3(flt.xc(2436,1),flt.xc(2436,2),flt.xc(2436,3),'o','MarkerSize',10,'Color','b');
    title('Dip-slip traction component');
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    set(gca,'view',[-227 42]);
    set(gca,'clim',[-1 1]*max((get(gca,'clim')))/1e0);
    axis equal
    set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3);
    
    subplot(2,1,2);cla;hold on;
    
    %toplot=evl.evt{1}.KL{1,5}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,5}*(evt{1}.slip.*sind(evt{1}.rake));
    %toplot=evl.evt{1}.KL{1,6}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,6}*(evt{1}.slip.*sind(evt{1}.rake));
    toplot=evl.evt{1}.KL{1,4}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,4}*(evt{1}.slip.*sind(evt{1}.rake));
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
        
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    title('Stress component s_{23}');
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    set(gca,'view',[-227 42]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
    axis equal
    set(gca,'xlim',[-80 350]*1e3,'ylim',[-150 100]*1e3);
    fprintf('Stress interaction: review before simulation\n');
    return
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% rheology
evl.shz.n=0*shz.n+1;
evl.shz.m=0*shz.n+1;
evl.shz.etaM=shz.etaM*0+1e19/(365*24*60*60)/1e6;
evl.shz.etaK=1e-1*evl.shz.etaM;
%evl.shz.etaM=1e0*arr(evl.shz.T0);
%evl.shz.etaK=1e-2*evl.shz.etaM;

evl.flt.Vo=0*evl.flt.l+3e0;
evl.flt.l=0*flt.l+7.5e-2;
evl.flt.Vpl=0*flt.l+0.001;
evl.flt.b=evl.flt.a-8e-2;

s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-20e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)-20e3)*sind(strike);

% pin patches with negative stress change at t=0
evl.flt.pinnedPosition=find((evl.evt{1}.KK{1,2}*(evl.evt{1}.src.slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evt{1}.rake)))<0);
evl.flt.pinnedPosition=sort(unique([flt.pinnedPosition; ...
    find(d<0)]));

% initial conditions
y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);

% integration options
%options=ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','on','oDir','./postseismic_pinned');
options=unicycle.ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','off');
tic
% provides evl.t and evl.y
tpts=unique(sort([gpsdata{1}(:,1); gpsdata{2}(:,1); gpsdata{3}(:,1); gpsdata{4}(:,1)]));
%tpts=[0:0.0027:0.7]
evl.ode45(tpts,y0,options);
toc

%% surface data simulation

% plot GPS time series
if true
    colors='mbrkcy';
    % GPS data (network.dat contains a list of stations with name and coordinates)
    if ~exist('gps','var')
        tic
        gps=unicycle.manifold.gpsReceiver('./gps/caltech-dase-network.dat',evl,3);
        toc
    else
        gps.simulation(evl);
    end
    ylim=[-0.1 0.1];
    figure(1);clf;set(gcf,'Name','GPS time series')
    index={1,2,3,4};
    for i=1:12
    subplot(4,3,i);cla;
    end
    for k=1:length(index)
        stationId=index{k};
        subplot(4,3,1+3*(k-1));
        hold on, box on;
        toplota=gps.ua((stationId-1)*gps.vecsize+2,:)'*cosd(strike)+...
            gps.ua((stationId-1)*gps.vecsize+1,:)'*sind(strike);
        toplotv=gps.uv((stationId-1)*gps.vecsize+2,:)'*cosd(strike)+...
            gps.uv((stationId-1)*gps.vecsize+1,:)'*sind(strike);
        %toplota=gps.ua((stationId-1)*gps.vecsize+1,:)';
        %toplotv=gps.uv((stationId-1)*gps.vecsize+1,:)';
        plot(evl.t,toplota,'-.','lineWidth',1,'color',colors(k));
        plot(evl.t,toplotv,'--','lineWidth',1,'color',colors(k));
        plot(evl.t,toplota+toplotv,'b-','lineWidth',2,'color',colors(k));
        toplotd=gpsdata{k}(:,3)*cosd(strike)+...
            gpsdata{k}(:,2)*sind(strike);
        plot(gpsdata{k}(:,1),toplotd.*1e-3,'b-','lineWidth',5,'color',colors(k));
        title(sprintf('station %s, Fault-parallel displacement',gps.stationName{stationId}));
        ylabel('Fault-parallel displacement (m)');
        set(gca,'ylim',ylim)
        
        subplot(4,3,2+3*(k-1));
        hold on, box on;
        toplota=-(gps.ua((stationId-1)*gps.vecsize+2,:)'*cosd(azimuth)+...
            gps.ua((stationId-1)*gps.vecsize+1,:)'*sind(azimuth));
        toplotv=-(gps.uv((stationId-1)*gps.vecsize+2,:)'*cosd(azimuth)+...
            gps.uv((stationId-1)*gps.vecsize+1,:)'*sind(azimuth));
        %toplota=gps.ua((stationId-1)*gps.vecsize+2,:)';
        %toplotv=gps.uv((stationId-1)*gps.vecsize+2,:)';
        plot(evl.t,toplota,'-.','lineWidth',1,'color',colors(k));
        plot(evl.t,toplotv,'--','lineWidth',1,'color',colors(k));
        plot(evl.t,toplota+toplotv,'b-','lineWidth',2,'color',colors(k));
        toplotd=-(gpsdata{k}(:,3)*cosd(azimuth)+...
            gpsdata{k}(:,2)*sind(azimuth));
        plot(gpsdata{k}(:,1),toplotd.*1e-3,'b-','lineWidth',5,'color',colors(k));
        title(sprintf('Station %s, Fault-perpendicular displacement',gps.stationName{stationId}));
        ylabel('Fault-perpendicular displacement (m)');
        set(gca,'ylim',ylim)
        
        subplot(4,3,3+3*(k-1));
        hold on, box on;
        toplota=gps.ua((stationId-1)*gps.vecsize+3,:)';
        toplotv=gps.uv((stationId-1)*gps.vecsize+3,:)';
        plot(evl.t,toplota,'-.','lineWidth',1,'color',colors(k));
        plot(evl.t,toplotv,'--','lineWidth',1,'color',colors(k));
        plot(evl.t,toplota+toplotv,'b-','lineWidth',2,'color',colors(k));
        plot(gpsdata{k}(:,1),gpsdata{k}(:,4).*1e-3,'b-','lineWidth',5,'color',colors(k));
        title(sprintf('Station %s, Vertical displacement',gps.stationName{stationId}));
        ylabel('Vertical displacement (m)'), xlabel('time (yr)');
        set(gca,'ylim',ylim)
    end

end


%% map

% plot GPS in map view
if true
    colormap default
    scale=2e6;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
   
    timeIndex=ceil(length(evl.t));
    %shz.plotShearZoneIndex();
    shz.plotShearZoneEdges();
    plot(gps.x(:,1),gps.x(:,2),'^');
    gpsdataplot=zeros(length(gpsdata),5);
    for i=1:length(gpsdata)
        stationId=index{i};
        timedataIndex=dsearchn(gpsdata{i}(:,1),evl.t(timeIndex))
        gpsdataplot(i,:)=[gps.x(stationId,1) gps.x(stationId,2) ...
            1e-3*gpsdata{i}(timedataIndex,2) 1e-3*gpsdata{i}(timedataIndex,3) ...
            1e-3*gpsdata{i}(timedataIndex,4)];
        1e-3*gpsdata{i}(timedataIndex,4)
    end
    scatter(gpsdataplot(:,1),gpsdataplot(:,2),700,gpsdataplot(:,5),'filled')
     quiver(gpsdataplot(:,1),gpsdataplot(:,2), ...
        scale*gpsdataplot(:,3),scale*gpsdataplot(:,4),0,'r','lineWidth',3);
    
   % scale*gps.ua(1:gps.vecsize:end,timeIndex)
    %evt{1}.plotPatch()
    scatter(gps.x(:,1),gps.x(:,2),250,gps.ua(3:gps.vecsize:end,timeIndex),'filled')
    quiver(gps.x(:,1),gps.x(:,2), ...
           scale*(gps.ua(1:gps.vecsize:end,timeIndex)+gps.uv(1:gps.vecsize:end,timeIndex)), ...
           scale*(gps.ua(2:gps.vecsize:end,timeIndex)+gps.uv(2:gps.vecsize:end,timeIndex)),0,'b','lineWidth',2);
       

       
       % add in quiver command for data vectors as well as model
    
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
    end
    
    stressComponent=3;
    %toplot=evl.evt{1}.KL{1,stressComponent}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,stressComponent}*(evt{1}.slip.*sind(evt{1}.rake));
    %shz.plotShearZoneEdges(toplot);
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight

    set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3);
    h=colorbar();
    ylabel(h,'Uplift (m)')
    title(sprintf('Surface displacements at time %f yr',evl.t(timeIndex)));
end

%% static properties 3d plot

if true
        %figure(1001);clf;set(gcf,'name','Stress interactions');
    %subplot(2,1,1);gca;hold on
    figure(8);clf;set(gcf,'name','Postseismic stress interactions')
       subplot(1,2,1);gca;hold on 
    %set(gcf,'name','snapshot')
    
    timeIndex=ceil(length(evl.t)/1);
    
    %evt{1}.plotPatch();
    
    % strain in shear zone
  %  toplot=evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,timeIndex)*1e6;
    
    % stress in shear zone
 %   toplot=evl.y(evl.flt.N*evl.flt.dgf+7:evl.shz.dgf:end,timeIndex);
    
   % shz.plotShearZoneEdges(toplot);
   % shz.plotShearZoneEdges();
    
    % afterslip
    slip1=evl.evt{1}.src.slip;
    toplot=evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex)%+slip1;
    %toplot=evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex);%+slip1;
    
    flt.plotPatch(toplot);
    flt.plotPatch();
    
	colormap(jet)
    h=colorbar('SouthOutside');
    ylabel(h,'Slip (m)')
    
    title(sprintf('Afterslip (m) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;view([94 46]);%view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[0 1]*max(abs(clim(:))));
    %set(gca,'clim',[-1 1]*5e-2);
    %box on, view([-227 42]), xlabel('x (km)'), ylabel('y (km)')
    grid on
    
    
    
    subplot(1,2,2); cla; hold on
      timeIndex=50;
        set(gcf,'name','snapshot')
        toplot=evl.y(evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end,timeIndex)*1e6;
        shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    	colormap(jet)
    h=colorbar('SouthOutside');
    ylabel(h,'Strain')
        title(sprintf('Strain component e_{23} t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;view([94 46]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:)))/1e1);
end

%% slip evolution movie

if true
    isVideoSave=false;
    if isVideoSave
        vidObj = VideoWriter('viscoelastic_1/strain_1.avi');
        open(vidObj);
    end
    % bounds
    xlim=[-80 220]*1e3;
    ylim=[-150 100]*1e3;
    component=2; % dip slip
    
    figure(9);clf;
    view([-227 42]);
    for k=2:25:length(evl.t)
        
        hold on
        
        % strain component e12 (microstrain)
        toplot=evl.y((evl.flt.N*evl.flt.dgf)+5:evl.shz.dgf:end,k)*1e6;
        
        % stress component e13 (MPa)
        %toplot=evl.y(11:evl.shz.dgf:end,k);
        
        shz.plotShearZoneEdges(toplot);
        shz.plotShearZoneEdges();
        %evt{1}.plotPatch();
        
        % afterslip in dip direction (m)
        toplot=evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k);
        flt.plotPatch(toplot);
        flt.plotPatch();
        
        plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2);
        
        colormap(gca,parula)
        h=colorbar();
        ylabel(h,'strain')
        
        title(sprintf('t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight; box on, grid on;
        
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-65 0]*1e3);
        set(gca,'clim',[-1 1]*6);
        box on
        
        xlabel('x (km)'), ylabel('y (km)');
        %set(gca,'nextplot','replacechildren');
        drawnow
        if isVideoSave
            currentFrame=getframe; 
            writeVideo(vidObj,currentFrame);
        end
    end
    if isVideoSave
        close(vidObj);
    end
end

%% Exporting

timeIndex=ceil(length(evl.t)/1);
e11=evl.y(evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end,timeIndex)*1e6;
e12=evl.y(evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end,timeIndex)*1e6;
e13=evl.y(evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end,timeIndex)*1e6;
e22=evl.y(evl.flt.N*evl.flt.dgf+4:evl.shz.dgf:end,timeIndex)*1e6;
e23=evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,timeIndex)*1e6;
e33=evl.y(evl.flt.N*evl.flt.dgf+6:evl.shz.dgf:end,timeIndex)*1e6;
evl.shz.exportSHZ('test_output2.shz',e11,e12,e13,e22,e23,e33);
