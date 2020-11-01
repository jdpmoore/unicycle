% Script POSTSEISMIC_1 simulates postseismic deformation following the 
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
% The equation of state for the flow of quartzites at high temperature is
% used to predict the viscoelastic relaxation of the lower crust. The
% geometry is described in the file
%
%   ./faults/lowercrust.shz
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
% AUTHOR: 
% Sylvain Barbot (July 14, 2016), Earth Observatory of Singapore

%clear all
%close all

addpath /Users/James/Documents/src/unicycle/matlab
import unicycle.*

load './gmt/boundaries.dat';

s2y=60*60*24*365;
y2s=1./s2y;
colors='kmcrgby';
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

%flt=geometry.triangleReceiver('./faults/qiu+15_1_receiver',greens.nikkhoo15(G,nu));
%flt=geometry.triangleReceiver('./faults/test',greens.nikkhoo15(G,nu));
%flt=geometry.receiver('./faults/receiver',greens.okada92(G,nu));
flt=geometry.receiver('./faults/feng+15_1_receiver',greens.okada92(G,nu));


% rate-and-state coefficients
flt.a=0*flt.l+1e-2;
flt.b=flt.a-4e-3;
strike=-70;
s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-0e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)+0e3)*sind(strike);
%index=find(1==omega(s/120e3).*omega(d/140e3));
index=find(1==omega(s/120e3).*omega(flt.xc(:,3)/40e3));
%flt.b(index)=flt.a(index)+4e-3;

% reference velocity (m/yr)
flt.Vo=0*flt.l+1e+1;

% characteristic slip distance
flt.l=0*flt.l+7.5e-2;

% confining pressure (MPa)
flt.sigma=0*flt.l+500;

% plate loading
flt.Vpl=0*flt.l+0.02;
azimuth=-160;
azimuthVector=[0*flt.l+sind(azimuth),0*flt.l+cosd(azimuth),0*flt.l];
flt.Vrake=atan2(sum(flt.dv.*azimuthVector,2),sum(flt.sv.*azimuthVector,2))*180/pi;
flt.rake=flt.Vrake;

% patches of interest
% flt.observationPoints={ ...
%     unicycle.geometry.observationPoint(0902), ... % trench
%     unicycle.geometry.observationPoint(2070), ... % Gorkha décollement
%     unicycle.geometry.observationPoint(4688), ... % upper décollement
%     unicycle.geometry.observationPoint(6564) ...  % far field
%     }; % for MFT_receiver_v20_short

%flt.observationPoints={ ...
%    unicycle.geometry.observationPoint(3) % trench
%    }; % for MFT_receiver_v20_short

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=[];%geometry.source();

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%evt={geometry.coseismicTriangle('./faults/qiu+15_1',0,greens.nikkhoo15(G,nu))};
evt={geometry.coseismicPatch('faults/feng+15_1_patch.flt',0,greens.okada92(G,nu))};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%      R E C E I V E R   S H E A R   Z O N E S         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

shz=geometry.shearZoneReceiver('./faults/lowercrust',greens.shearZone16(G,nu));

s=+(shz.x(:,2)+250e3)*cosd(strike)+(shz.x(:,1)-0e3)*sind(strike);
d=+(shz.x(:,1)+140e3)*cosd(strike)-(shz.x(:,2)+0e3)*sind(strike);

ramp=@(x) x.*omega(x-1/2)+heaviside(x-1);

shz.x(:,3)=shz.x(:,3)-35e3*ramp(d/1e5);
shz.xc(:,3)=shz.xc(:,3)-35e3*ramp(d/1e5);

% rheology
shz.n=0*shz.n+1;
shz.etaM=zeros(shz.N,1)+3.171e3; % MPa yr

% plate loading
shz.epsilonPlate=[0,1,0,0,0,0]*1e-15;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            G E O M E T R Y   C H E C K               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & properties');
    hold on
    
    shz.plotShearZoneEdges();
    %shz.plotUnitVectors(1e3);
    
    sc=1e4;
    %flt.plotPatch(flt.b-flt.a);
    flt.plotPatch();
    %flt.plotPatchIndex();
    evt{1}.plotPatch(evt{1}.slip);
    
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100)
    end
    
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    colormap(jet)
    h=colorbar();
    ylabel(h,'friction properties');
    box on, grid on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('(b-a) frictional properties')
    %set(gca,'view',[-227 42]);
    axis equal
    set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3)
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
    evl=ode.rateStrengtheningMaxwell([],flt,shz,evt);
    toc
else
    % skip building stress kernel and update velocities & friction properties
    evl.flt=flt;
    evl.flt.dgf=4;
    evl.shz=shz;
    evl.shz.dgf=12;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%          S T R E S S   I N T E R A C T I O N         %
%               V I S U A L I Z A T I O N              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1001);clf;set(gcf,'name','Stress interactions');
    subplot(2,1,1);gca;hold on
    
    %shz.plotShearZoneEdges();
    
    toplot=evl.evt{1}.KK{1,2}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(evt{1}.slip.*sind(evt{1}.rake));
    flt.plotPatch(toplot);
    
    sc=1e4;
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    set(gca,'view',[-227 42]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/5e0);
    axis equal
    set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3);
    fprintf('Stress interaction: review before simulation\n');
    
    subplot(2,1,2);cla;hold on;
    
    %toplot=evl.evt{1}.KL{1,5}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,5}*(evt{1}.slip.*sind(evt{1}.rake));
    toplot=evl.evt{1}.KL{1,6}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,6}*(evt{1}.slip.*sind(evt{1}.rake));
    shz.plotShearZoneEdges(toplot);
    
    %flt.plotPatch();
    
    shz.plotShearZoneEdges();
    
    sc=1e4;
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
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
evl.flt.isRakeConstraint=true;
% initial conditions
y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-13,'InitialStep',1e-11,'Stats','on','oDir','./postseismic_2');
tic
% provides evl.t and evl.y
evl.ode45([0 3],y0,options);
toc

%% export to Paraview

if false
    evl.exportvtp(evl.t,evl.y,1e0,50,'./run4/vtk');
end

%% surface data simulation

% plot GPS time series
if true
    % GPS data (network.dat contains a list of stations with name and coordinates)
    if ~exist('gps','var')
        gps=unicycle.manifold.gpsReceiver('./gps/caltech-dase-network.dat',evl,3);
    else
        gps.simulation(evl);
    end
    
    figure(1);clf;set(gcf,'Name','GPS time series')
    index={4,5};
    subplot(3,1,1);cla;
    subplot(3,1,2);cla;
    subplot(3,1,3);cla;
    
    for k=1:length(index)
        stationId=index{k};
        subplot(3,1,1);
        hold on, box on;
        toplot=gps.uv((stationId-1)*gps.vecsize+2,:)'*cosd(strike)+...
            gps.uv((stationId-1)*gps.vecsize+1,:)'*sind(strike);
        %toplot=gps.uv((stationId-1)*gps.vecsize+1,:)';
        plot(evl.t,toplot,'b-','lineWidth',1,'color',colors(stationId));
        title(sprintf('station %s, Fault-parallel displacement',gps.stationName{stationId}));
        ylabel('East displacement (m)');
        
        subplot(3,1,2);
        hold on, box on;
        toplot=-(gps.uv((stationId-1)*gps.vecsize+2,:)'*cosd(azimuth)+...
            gps.uv((stationId-1)*gps.vecsize+1,:)'*sind(azimuth));
        %toplot=gps.uv((stationId-1)*gps.vecsize+2,:)'
        plot(evl.t,toplot,'b-.','lineWidth',1,'color',colors(stationId));
        title(sprintf('Station %s, Fault-perpendicular displacement',gps.stationName{stationId}));
        ylabel('North displacement (m)');
        
        subplot(3,1,3);
        hold on, box on;
        toplot=gps.uv((stationId-1)*gps.vecsize+3,:)';
        plot(evl.t,toplot,'b-','lineWidth',1,'color',colors(stationId));
        title(sprintf('Station %s, Vertical displacement',gps.stationName{stationId}));
        ylabel('Vertical displacement (m)'), xlabel('time (yr)');
    end
end


%% map

% plot GPS in map view
if true
    scale=2e6;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    
    timeIndex=length(evl.t);
    %shz.plotShearZoneIndex();
    shz.plotShearZoneEdges();
    plot(gps.x(:,1),gps.x(:,2),'^');
    %evt{1}.plotPatch()
    scatter(gps.x(:,1),gps.x(:,2),250,gps.uv(3:gps.vecsize:end,timeIndex),'filled')
    quiver(gps.x(:,1),gps.x(:,2), ...
           scale*gps.uv(1:gps.vecsize:end,timeIndex), ...
           scale*gps.uv(2:gps.vecsize:end,timeIndex),0);
    
    for k=1:length(gps.stationName)
        text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
    end
    
    stressComponent=3;
    toplot=evl.evt{1}.KL{1,stressComponent}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,stressComponent}*(evt{1}.slip.*sind(evt{1}.rake));
    shz.plotShearZoneEdges(toplot);
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight

    %set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3);
    h=colorbar();
    ylabel(h,'Uplift (m)')
    title(sprintf('Surface displacements at time %f yr',evl.t(timeIndex)));
end

%% static properties 3d plot

if false
    figure(8);clf;
    
    set(gcf,'name','source and receiver fault geometry')
    hold on
    %timeIndex=ceil(length(evl.t)/2);
    [~,timeIndex1]=min((evl.t-100).^2);
    [~,timeIndex2]=min((evl.t-101).^2);
    Dt=evl.t(timeIndex2)-evl.t(timeIndex1);
    
    thrust_velovity=evl.y(1:evl.dgf:end,timeIndex2).*cosd(flt.Vrake)+...
                    evl.y(2:evl.dgf:end,timeIndex2).*sind(flt.Vrake)-...
                    evl.y(1:evl.dgf:end,timeIndex1).*cosd(flt.Vrake)-...
                    evl.y(2:evl.dgf:end,timeIndex1).*sind(flt.Vrake);
    
	% instantaneous strike- and dip-direction velocity
	vss=(evl.y(1:evl.dgf:end,timeIndex2)-evl.y(1:evl.dgf:end,timeIndex1))/Dt;
    vds=(evl.y(2:evl.dgf:end,timeIndex2)-evl.y(2:evl.dgf:end,timeIndex1))/Dt;
       
                
	% plunge in the cross-section parallel to Vrake
    plunge=atan2(cosd(flt.Vrake).*flt.sv(:,3)+sind(flt.Vrake).*flt.dv(:,3),sum((repmat(cosd(flt.Vrake),1,3).*flt.sv+repmat(sind(flt.Vrake),1,3).*flt.dv).*azimuthVector,2))*180/pi;
    
    % slip rate with unit horizontal component
    hors=cosd(flt.Vrake).*(evl.Kss*(vss)+evl.Kds*(vds))+ ...
         sind(flt.Vrake).*(evl.Ksd*(vss)+evl.Kdd*(vds));
 
    horn=evl.Ksn*(vss)+evl.Kdn*(vds);
    flt.plotPatch(hors+1*horn);
       
    %flt.plotPatch(dip);
      
    %flt.plotSlipVectors(evl.y(1:evl.dgf:end,timeIndex2)'-...
    %                    evl.y(1:evl.dgf:end,timeIndex1)',...
    %                    evl.y(2:evl.dgf:end,timeIndex2)'-...
    %                    evl.y(2:evl.dgf:end,timeIndex1)',2e3)
	colormap(jet)
    h=colorbar();
    ylabel(h,'slip (m)')
    for k=1:length(index)
        %plot3(flt.xc(index{k},1),flt.xc(index{k},2),flt.xc(index{k},3),'r+','MarkerSize',100)
    end
    
    title(sprintf('Coulomb stress rate (MPa/yr) t=%2.1e/%2.1e yr',evl.t(timeIndex1),evl.t(end)));
    axis equal tight;%view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    %set(gca,'clim',[-1 1]*max(abs(clim(:))));
    set(gca,'clim',[-1 1]*5e-2);
    box on, view([-227 42]), xlabel('x (km)'), ylabel('y (km)')
    grid on
end


%% static properties 3d plot

if true
    figure(8);clf;
    
    set(gcf,'name','snapshot')
    hold on
    
    timeIndex=ceil(length(evl.t)/1);
   
    %evt{1}.plotPatch();
    
    % strain
    toplot=evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,timeIndex)*1e6;
    
    % stress
    %toplot=evl.y(7:evl.shz.dgf:end,timeIndex);
    shz.plotShearZoneEdges(toplot/1e1);
    shz.plotShearZoneEdges();
    toplot=evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,timeIndex);
    flt.plotPatch(toplot);
    
	colormap(jet)
    h=colorbar();
    ylabel(h,'Stress (m)')
    
    title(sprintf('Stress component s_{12} (MPa) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;%view([-30,20]);
    %set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-145 0]);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    %set(gca,'clim',[-1 1]*5e-2);
    box on, view([-227 42]), xlabel('x (km)'), ylabel('y (km)')
    grid on
end

%% slip evolution movie

if false
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
    for k=2:50:length(evl.t)
        
        hold on
        
        % strain component e12 (microstrain)
        %toplot=evl.y(2:evl.shz.dgf:end,k)*1e6;
        
        % stress component e13 (MPa)
        toplot=evl.y(evl.flt.N*evl.flt.dgf+11:evl.shz.dgf:end,k);
        
        shz.plotShearZoneEdges(toplot);
        shz.plotShearZoneEdges();
        evt{1}.plotPatch();
        
        plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2);
        
        colormap(gca,parula)
        h=colorbar();
        ylabel(h,'strain')
        
        title(sprintf('t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight; box on, grid on;
        
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-65 0]*1e3);
        set(gca,'clim',[-1 1]*1);
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



