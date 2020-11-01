% Script POSTSEISMIC_2 simulates postseismic deformation following slip on 
% a strike slip fault. 
%
% AUTHOR: 
% Sylvain Barbot (July 16, 2016), Earth Observatory of Singapore

clear all
import unicycle.*

s2y=60*60*24*365;
y2s=1./s2y;
colors='kmcrgby';
G=30e3;
nu=1/4;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% homogeneous elastic half space, Poisson's solid
earthModel=greens.okada92(30e3,1/4);

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

index=find(1==omega((flt.xc(:,2)-100e3)/100e3).*omega((flt.xc(:,3)+20e3)/20e3));
flt.b(index)=flt.a(index)+4e-3;

% reference velocity (m/yr)
flt.Vo=0*flt.L+1e+1;

% characteristic slip distance
flt.l=0*flt.L+0.3;

% confining pressure (MPa)
flt.sigma=0*flt.L+100;
flt.dgf=5;

fprintf('h* = %2.2e km\n',max(earthModel.G*flt.l./(flt.b-flt.a)./flt.sigma/1e3))

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

evt={geometry.coseismicPatch('./faults/strikeslip.flt',0,greens.okada92(G,nu))};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%      R E C E I V E R   S H E A R   Z O N E S         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

shz=geometry.shearZoneReceiver('./faults/lowercrust',greens.shearZone16(G,nu));

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                S A N I T Y   C H E C K               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','geometry & properties');
    hold on
    
    sc=1e4;
    shz.plotShearZoneEdges();
    %shz.plotUnitVectors(1e3);
    
    evt{1}.plotPatch();
    
    for k=1:length(shz.observationPoints)
        plot3( ...
            shz.observationPoints{k}.xc(1), ...
            shz.observationPoints{k}.xc(2), ...
            shz.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100)
    end
    
    colormap(jet)
    h=colorbar();
    ylabel(h,'properties');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Properties')
    %set(gca,'view',[-227 42]);
    axis equal
    set(gca,'xlim',[-80 80]*1e3,'ylim',[-80 80]*1e3)
    fprintf('Sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S T R E S S   K E R N E L S             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl')
    % build stress kernels for integral equation
    tic
    evl=ode.maxwell(evt,flt,shz,src);
    toc
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.LL{1,1},2)==size(evl.shz.xc,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.shz=shz;
    evl.shz.dgf=12;
    evl.evt=evt;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%          S T R E S S   I N T E R A C T I O N         %
%               V I S U A L I Z A T I O N              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1000);clf;set(gcf,'name','Stress interactions');
    hold on
    
    evt{1}.plotPatch();
    
    sc=1e4;
    
    % pressure interactions
    %toplot=evl.LL{1,1}+evl.LL{1,4}+evl.LL{1,6}+ ...
    %       evl.LL{4,1}+evl.LL{4,4}+evl.LL{4,6}+ ...
    %       evl.LL{6,1}+evl.LL{6,4}+evl.LL{6,6};
    
    % coseismic stress change
    toplot=evl.KL{1,2}*(evl.src.slip.*cosd(evl.src.rake));
    
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    %shz.plotUnitVectors(1e3);
    
    for k=1:length(shz.observationPoints)
        plot3( ...
            shz.observationPoints{k}.xc(1), ...
            shz.observationPoints{k}.xc(2), ...
            shz.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    set(gca,'view',[-227 42]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
    axis equal
    set(gca,'xlim',[-1 1]*60e3,'ylim',[-1 1]*60e3);
    fprintf('Stress interaction: review before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% initial conditions (other than coseismic stress perturbation)
y0=zeros(1,evl.shz.N*evl.shz.dgf);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-13,'InitialStep',1e-5,'Stats','on');
tic
% provides evl.t and evl.y
evl.ode45([0 10],y0,options);
toc


%% export to Paraview

if false
    evl.exportvtp(evl.t,evl.y,1e0,50,'./run4/vtk');
end

%% surface data simulation

% plot GPS time series
if true
    % GPS data (network.dat contains a list of stations with name and coordinates)
    gps=unicycle.manifold.gpsReceiver('./gps/gps-network.dat',evl,3);
    
    figure(1);clf;set(gcf,'Name','GPS time series')
    stationId=4;
    subplot(3,1,1);cla;
    hold on, box on;
    toplot=gps.uv((stationId-1)*gps.vecsize+2,:)';
    plot(evl.t,toplot,'b-','lineWidth',5);
    title(sprintf('station %s, north component',gps.stationName{stationId}));
    ylabel('East displacement (m)');
    
    subplot(3,1,2);cla;
    toplot=gps.uv((stationId-1)*gps.vecsize+1,:)';
    plot(evl.t,toplot,'b-','lineWidth',5);
    title(sprintf('station %s, east component',gps.stationName{stationId}));
    ylabel('North displacement (m)');
    
    subplot(3,1,3);cla;
    toplot=gps.uv((stationId-1)*gps.vecsize+3,:)';
    plot(evl.t,toplot,'b-','lineWidth',5);
    title(sprintf('station %s, vertical displacement',gps.stationName{stationId}));
    ylabel('Vertical displacement (m)'), xlabel('time (yr)');
end


%% map

% plot GPS in map view
if true
    scale=1e5;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    
    timeIndex=length(evl.t);
    %shz.plotShearZoneIndex();
    shz.plotShearZoneEdges();
    plot(gps.x(:,1),gps.x(:,2),'^');
    evt{1}.plotPatch()
    scatter(gps.x(:,1),gps.x(:,2),500,gps.uv(3:gps.vecsize:end,timeIndex),'filled')
    quiver(gps.x(:,1),gps.x(:,2), ...
           scale*gps.uv(1:gps.vecsize:end,timeIndex), ...
           scale*gps.uv(2:gps.vecsize:end,timeIndex),0);
    
    for k=1:length(gps.stationName)
        %text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
    end
    
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
    set(gca,'xlim',[-1 1]*60e3)
    h=colorbar();
    ylabel(h,'Uplift (m)')
    title(sprintf('Surface displacements at time %f yr',evl.t(timeIndex)));
end

%% static properties 3d plot

if true
    figure(8);clf;
    
    set(gcf,'name','snapshot')
    hold on
    
    timeIndex=3;%ceil(length(evl.t)/3);
    
    evt{1}.plotPatch();
    
    % strain
    %toplot=evl.y(6:evl.shz.dgf:end,timeIndex)*1e6;
    
    % stress
    toplot=evl.y(8:evl.shz.dgf:end,timeIndex);
    shz.plotShearZoneEdges(toplot);
    
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
    xlim=[-60 60]*1e3;
    ylim=[-60 60]*1e3;
    component=2; % dip slip
    
    figure(9);clf;
    view([-227 42]);
    for k=2:1:length(evl.t)
        
        hold on
        
        % strain component e12 (microstrain)
        %toplot=evl.y(2:evl.shz.dgf:end,k)*1e6;
        
        % stress component e12 (MPa)
        toplot=evl.y(9:evl.shz.dgf:end,k);
        
        shz.plotShearZoneEdges(toplot);
        shz.plotShearZoneEdges();
        evt{1}.plotPatch();
        
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



