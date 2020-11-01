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
% Sylvain Barbot and James Moore (July 29, 2016), Earth Observatory of Singapore

%clear all
%close all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

load './matlab/japan_coasts_km.dat';
load './matlab/japan_faults_km.dat';
load './matlab/faults_maps2.dat';
load './matlab/plate_boundaries_km.dat';
load './matlab/ryukyu_km.dat';
load './matlab/volcano_km.dat';
load './gmt/post_seismicity_km.dat'
load './Inversion/GMatrices/kernels.mat'
%formatSpec=;
gpsdata=readtable('./gps/disp_pos3m_gps.dat','Format','%d%s%f%f%f%f%f%f%f%f','ReadVariableNames',true);
[xm,ym,zm] = unicycle.export.grdread('./gmt/litho1.0_moho_km.grd');



%gpsdata=cell(4,1);
%gpsdata{1} = importdata('./gps/KKN4out.dat');
%gpsdata{2} = importdata('./gps/NASTout.dat');
%gpsdata{3} = importdata('./gps/CHLMout.dat');
%gpsdata{4} = importdata('./gps/SNDLout.dat');

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

%flt=geometry.receiver('./faults/kumamoto-hayes16-prel',greens.okada92(G,nu));
flt=geometry.receiver('./faults/swei16',greens.okada92(G,nu));

% rate-and-state coefficients
flt.a=0*flt.l+1e-2;
flt.b=flt.a-8e-3;

strike=mean(flt.strike);
s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-0e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)+0e3)*sind(strike);
%index=find(1==omega(s/120e3).*omega(d/140e3));
index=find(1==omega(s/120e3).*omega(flt.xc(:,3)/40e3));
%flt.b(index)=flt.a(index)+4e-3;

% reference velocity (m/yr)
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

%evt={geometry.coseismicPatch('./faults/kumamoto-hayes16-prel.flt',0,greens.okada92(G,nu))};
evt={geometry.coseismicPatch('./faults/swei16.flt',0,greens.okada92(G,nu))};


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%      R E C E I V E R   S H E A R   Z O N E S         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

shz=geometry.shearZoneReceiver('./faults/lc+sz+mw',greens.shearZone16(G,nu));



pos=1:shz.levels{1}.nShearZone;
strike=270;
dz=interp2(1e3*xm,1e3*ym,-1e3*zm,shz.x(pos,1),shz.x(pos,2));
s=+(shz.x(pos,2)+50e3)*cosd(strike)+(shz.x(pos,1)-0e3)*sind(strike);
d=+(shz.x(pos,1)-20e3)*cosd(strike)-(shz.x(pos,2)+0e3)*sind(strike);

ramp=@(x) x.*omega(x-1/2)+heaviside(x-1);
shz.x(pos,3)=shz.x(pos,3)+dz+shz.levels{1}.W;
shz.xc(pos,3)=shz.xc(pos,3)+dz+shz.levels{1}.W;

%shz.x(pos,3)=shz.x(pos,3)-50e3*ramp(d/1e5);
%shz.xc(pos,3)=shz.xc(pos,3)-50e3*ramp(d/1e5);

% rheology
shz.n=0*shz.n+3;
shz.m=0*shz.n+2.7;
%shz.etaM=shz.etaM*0+1e20/(365*24*60*60)/1e6;
%shz.T0 = importdata('temps.dat');
shz.etaM=shz.etaM*0+1e19/(365*24*60*60)/1e6;%arr(shz.T0);
%shz.etaM=shz.etaM.*1e1.^(-ramp(d/1e5));
shz.etaK=1e-3*shz.etaM; %.*1e10;
%shz.etaK=shz.etaM*0+1e16/(365*24*60*60)/1e6;


% plate loading
shz.epsilonPlate=[0,1,0,0,0,0]*1e-15;
%out=shz.x(:,3), save('./testdepths.txt','out','-ascii')
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            G E O M E T R Y   C H E C K               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if true
    figure(1000);clf;set(gcf,'name','geometry & properties');
    hold on
  
         %shz.plotShearZoneEdges(log10(3.1536e+13*shz.etaM));
    %    shz.plotShearZoneEdges();
    
    %shz.plotUnitVectors(1e3);
    %slip1=evl.evt{1}.src.slip;slip1(slip1<0.4)=NaN;
    %flt.plotPatch(slip1);
    
    %sc=1e4;
    %flt.plotPatch(flt.b-flt.a);
    
    %flt.plotPatchIndex();
    evt{1}.plotPatch(evt{1}.slip);
    
      % shz.plotShearZoneEdges(log10(3.1536e+13*shz.etaM));
      shz.plotShearZoneEdges();
       % shz.plotShearZoneById((shz.levels{2}.starti+1):shz.levels{3}.starti);
    %shz.plotUnitVectors(1e3);
    
    sc=1e4;
   % flt.plotPatch(flt.b-flt.a);
    flt.plotPatch();
   % flt.plotPatchIndex();
 
    %evt{1}.plotPatch(evt{1}.slip);
    %evt{1}.plotSlipVectors(evt{1}.slip.*cosd(evt{1}.rake),evt{1}.slip.*sind(evt{1}.rake),3e3);
 %   flt.plotPatch();
    
  %  for k=1:length(flt.observationPoints)
  %      plot3( ...
  %          flt.observationPoints{k}.xc(1), ...
  %          flt.observationPoints{k}.xc(2), ...
  %          flt.observationPoints{k}.xc(3), ...
  %          [colors(k) '+'],'MarkerSize',100)
  %  end
    
    plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',2)
    % plot(1e3*japan_faults_km(:,1),1e3*japan_faults_km(:,2),'r-','lineWidth',2)
     plot(1e3*faults_maps2(:,1),1e3*faults_maps2(:,2),'r-','lineWidth',2)
 %    plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
   %  plot(syngpsx,syngpsy,'k^','markerfacecolor','b')
     plot3(1e3*ryukyu_km(:,1),1e3*ryukyu_km(:,2),1e3*ryukyu_km(:,3),'b-','lineWidth',2)
     % plot3(1e3*post_seismicity_km(:,1),1e3*post_seismicity_km(:,2),-1e3*post_seismicity_km(:,3),'ko','markerfacecolor','g')
  %   surf(1e3*xm,1e3*ym,-1e3*zm)
    colormap jet
    h=colorbar('southoutside');
   
    ylabel(h,'Coseismic slip (m)');
    box on, grid on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
  %  title('(b-a) frictional properties')
%    set(gca,'view',[-227 42]);
 %   set(gca,'view',[69.6000   -0.4000]);
    axis equal
    set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
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
    evl=ode.rateStrengtheningPower([],flt,shz,evt);
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

if true
      
    figure(1001);clf;set(gcf,'name','Stress interactions');
    subplot(1,2,1);gca;hold on
    slip=evl.evt{1}.src.slip;
    %shz.plotShearZoneEdges();
   toplot=-(evl.evt{1}.KK{1,1}*(slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,1}*(slip.*sind(evt{1}.rake)));
    %toplot(toplot>8)=8;
    %toplot=slip;
    
    %toplot=evl.KK{2,2};toplot=toplot(210,:);minmax(toplot)
    flt.plotPatch(toplot);
    %shz.plotShearZoneEdges(toplot);
    %flt.plotPatchIndex();
    %plot3(flt.xc(2436,1),flt.xc(2436,2),flt.xc(2436,3),'o','MarkerSize',10,'Color','b');
    title('Strike-slip traction component');
    sc=1e4;
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
         plot(1e3*japan_faults_km(:,1),1e3*japan_faults_km(:,2),'r-','lineWidth',2)
      plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
    %  plot(gps.x(:,1),gps.x(:,2),'kv','markerfacecolor','w')
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress (M Pa)');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    %set(gca,'view',[-227 42]);
    set(gca,'clim',[0 1]*max((get(gca,'clim')))/1e0);
     axis equal
    set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
    fprintf('Stress interaction: review before simulation\n');
        subplot(1,2,2);cla;hold on;
    
    %toplot=evl.evt{1}.KL{1,5}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,5}*(evt{1}.slip.*sind(evt{1}.rake));
    %toplot=evl.evt{1}.KL{1,6}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,6}*(evt{1}.slip.*sind(evt{1}.rake));
    
    export=cell(6,1);
    m=zeros(size(G_fault,2)+size(G_shrz,2)+size(G_poro,2)+size(G_ifault,2)+size(G_ishrz,2)+size(G_iporo,2),1);
    for k=1:6
        export{k}=evl.evt{1}.KL{1,k}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,k}*(evt{1}.slip.*sind(evt{1}.rake));
        m((1+size(G_fault,2)+(k-1)*shz{1}.N):6:(size(G_fault,2)+size(G_shrz,2)))=export{k}./1e-18;
    end
    poro=shz.levels{2}.starti:shz.levels{3}.starti;
    skk=(export{1}+export{4}+export{6})/2;
   export{1}(poro)=skk(poro);
   for k=2:6
       export{k}(poro)=NaN;
   end
    shz.exportSHZ('./faults/lowercrust_stress_swei16_poro.shz',export{:});
    J2=sqrt(((export{1}-export{4}).^2+(export{4}-export{6}).^2+(export{6}-export{1}).^2)./6+export{2}.^2+export{3}.^2+export{5}.^2);
    toplot=J2;
    %toplot=evl.evt{1}.KL{1,6}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,6}*(evt{1}.slip.*sind(evt{1}.rake));
    shz.plotShearZoneEdges(toplot);
    
    flt.plotPatch();
    
    shz.plotShearZoneEdges();
    
    sc=1e4;
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    
    plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
      %   plot(1e3*japan_faults_km(:,1),1e3*japan_faults_km(:,2),'g-','lineWidth',2)
      plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
    title('Stress component s_{22}');
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress (M Pa)');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    %set(gca,'view',[-227 42]);
    %set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/5e2);
        set(gca,'clim',[0 1]*max((get(gca,'clim')))/1e2);
    axis equal
set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
    fprintf('Stress interaction: review before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%        S A V E   S T R E S S    K E R N E L S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if false
    sKK=evl.KK;
    sKL=evl.KL;
    sLK=evl.LK;
    sLL=evl.LL;
    sFK=evl.FK;
    sFL=evl.FL;
    save('./StressKernels/KK.mat','sKK');
    save('./StressKernels/KL.mat','sKL');
    save('./StressKernels/LK.mat','sLK');
    save('./StressKernels/LL.mat','sLL');
    save('./StressKernels/FK.mat','sFK');
    save('./StressKernels/FL.mat','sFL');    
end

if false
   for i=1:length(evl.evt)
      fnameKK=['./StressKernels/evt_' num2str(i) '_KK.mat'];
      fnameKL=['./StressKernels/evt_' num2str(i) '_KL.mat'];
      sKK=evl.evt{i}.KK;
      sKL=evl.evt{i}.KL;
      save(fnameKK,'sKK');
      save(fnameKL,'sKL');
   end       
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% rheology
%evl.shz.T0 = importdata('temps2.dat');
evl.shz.n=0*shz.n+3;
%evl.shz.n=0*shz.n+1;
evl.shz.m=0*shz.n+2.7;
%evl.shz.m=0*shz.n+1;
evl.shz.etaM=shz.etaM*0+1e18/(365*24*60*60)/1e6;
%evl.shz.etaM(poro)=evl.shz.etaM(poro)*1e20;
%evl.shz.etaM=2e-5*arr(evl.shz.T0);
%evl.shz.etaK=shz.etaM*0+1e17/(365*24*60*60)/1e6;
%evl.shz.etaK=1e0*arr(evl.shz.T0);
evl.shz.etaK=5e-2*evl.shz.etaM;
evl.flt.Vo=0*evl.flt.l+3e1;
flt.l=0*flt.l+7.5e-2;

%s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-20e3)*sind(strike);
%d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)-20e3)*sind(strike);

% pin patches with negative stress change at t=0
flt.pinnedPosition=find((evl.evt{1}.KK{1,2}*(evl.evt{1}.src.slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evt{1}.rake)))<0);
%flt.pinnedPosition=sort(unique([flt.pinnedPosition; ...
%    find(d<0)]));

% initial conditions
y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);

% integration options
%options=ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','on','oDir','./postseismic_pinned');
options=ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','off');
tic
% provides evl.t and evl.y
%tpts=unique(sort([gpsdata{1}(:,1); gpsdata{2}(:,1); gpsdata{3}(:,1); gpsdata{4}(:,1)]));
day=0.0027;
tpts=[0:day:91*day];
evl.ode45(tpts,y0,options);
toc



%% export to Paraview

if false
    evl.exportvtp(evl.t,evl.y,1e0,50,'./run4/vtk');
end

%% surface data simulation


% plot GPS time series
if false
    colors='mbrkcy';
    % GPS data (network.dat contains a list of stations with name and coordinates)
    if ~exist('gps','var')
        tic
        gps=unicycle.manifold.gpsReceiver('./gps/kumamoto-gps.dat',evl,3);
        toc
    else
        gps.simulation(evl);
    end
    return
    ylim=[-0.1 0.1];
    figure(1);clf;set(gcf,'Name','GPS time series')
    index={4,5,13,15};
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


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              NORM OF DATA             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if false
    datL2=zeros(length(gpsdata),1);
    for i=1:length(gpsdata)
        stationId=index{i};
        [C,ia,ib] = intersect(gpsdata{i}(end,1),tpts);
        
        dat1=1e-3*(-(gpsdata{i}(end,3)*cosd(azimuth)+...
            gpsdata{i}(end,2)*sind(azimuth)))
        mod1=-(gps.ua((stationId-1)*gps.vecsize+2,ib)'*cosd(azimuth)+...
            gps.ua((stationId-1)*gps.vecsize+1,ib)'*sind(azimuth))...
            -(gps.uv((stationId-1)*gps.vecsize+2,ib)'*cosd(azimuth)+...
            gps.uv((stationId-1)*gps.vecsize+1,ib)'*sind(azimuth))
        datL2(i)=(dat1-mod1)^2;
        %modelN=gps.ua((stationId-1)*gps.vecsize+1,ib)+gps.uv((stationId-1)*gps.vecsize+1,ib);
        %modelE=gps.ua((stationId-1)*gps.vecsize+2,ib)+gps.uv((stationId-1)*gps.vecsize+2,ib);
        %modelD=gps.ua((stationId-1)*gps.vecsize+3,ib)+gps.uv((stationId-1)*gps.vecsize+3,ib);
        %datL2(i)=sqrt((modelN-gpsdata{i}(end,2).*1e-3).^2+(modelE-gpsdata{i}(end,3).*1e-3).^2+(modelD-gpsdata{i}(end,4).*1e-3).^2);
    end
    sqrt(sum(datL2))
end

%% map

% plot GPS in map view
if true
    if ~exist('gps','var')
        tic
        gps=unicycle.manifold.gpsReceiver('./gps/kumamoto-gps.dat',evl,3);
    %    gps=unicycle.manifold.gpsReceiver(syngps,evl,3);
        toc
    else
        gps.simulation(evl);
    end
    
            
    colormap default
    scale=2e6;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
   
    timeIndex=ceil(length(evl.t));
    %shz.plotShearZoneIndex();
    %shz.plotShearZoneEdges();
    plot(gps.x(:,1),gps.x(:,2),'^');
        plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
     %            plot(1e3*japan_faults_km(:,1),1e3*japan_faults_km(:,2),'r-','lineWidth',2)
     % plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
            axis equal

        
    
  %  gpsdataplot=zeros(length(gpsdata),5);
%     for i=1:length(gpsdata)
%         stationId=index{i};
%         timedataIndex=dsearchn(gpsdata{i}(:,1),evl.t(timeIndex))
%         gpsdataplot(i,:)=[gps.x(stationId,1) gps.x(stationId,2) ...
%             1e-3*gpsdata{i}(timedataIndex,2) 1e-3*gpsdata{i}(timedataIndex,3) ...
%             1e-3*gpsdata{i}(timedataIndex,4)];
%         1e-3*gpsdata{i}(timedataIndex,4)
%     end
     scatter(gpsdata.y,gpsdata.x,200,gpsdata.dU,'filled')
     
%     
   % scale*gps.ua(1:gps.vecsize:end,timeIndex)
    %evt{1}.plotPatch()
    scatter(gps.x(:,1),gps.x(:,2),50,gps.ua(3:gps.vecsize:end,timeIndex)+gps.uv(3:gps.vecsize:end,timeIndex),'filled')
     quiver(gpsdata.y,gpsdata.x, ...
         scale*gpsdata.dE,scale*gpsdata.dN,0,'r','lineWidth',1);
     quiver(gps.x(:,1),gps.x(:,2), ...
            scale*(gps.ua(1:gps.vecsize:end,timeIndex)+gps.uv(1:gps.vecsize:end,timeIndex)), ...
            scale*(gps.ua(2:gps.vecsize:end,timeIndex)+gps.uv(2:gps.vecsize:end,timeIndex)),0,'b','lineWidth',2);
        

       
       % add in quiver command for data vectors as well as model
    
%     for k=1:length(gps.stationName)
%         text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
%     end
    
    stressComponent=3;
    %toplot=evl.evt{1}.KL{1,stressComponent}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,stressComponent}*(evt{1}.slip.*sind(evt{1}.rake));
    %shz.plotShearZoneEdges(toplot);
    %plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight

   set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
    h=colorbar();
    ylabel(h,'Uplift (m)')
    title(sprintf('Surface displacements at time %d days',evl.t(timeIndex)/day));
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

if false
        %figure(1001);clf;set(gcf,'name','Stress interactions');
    %subplot(2,1,1);gca;hold on
    figure(8);clf;set(gcf,'name','Postseismic stress interactions')
     %  subplot(1,2,1);gca;hold on 
    %set(gcf,'name','snapshot')
    

    timeIndex=50;%ceil(length(evl.t)/1);
    
    %evt{1}.plotPatch();
    
    % strain in shear zone
  %  toplot=evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,timeIndex)*1e6;
    
    % stress in shear zone
 %   toplot=evl.y(evl.flt.N*evl.flt.dgf+7:evl.shz.dgf:end,timeIndex);
    
   % shz.plotShearZoneEdges(toplot);
   % shz.plotShearZoneEdges();
    
    % afterslip
    slip1=evl.evt{1}.src.slip;
    toplot=evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex)+slip1;
    
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
    
    return
    
    subplot(1,2,2); cla; hold on
      timeIndex=50;
        set(gcf,'name','snapshot')
        toplot=evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,timeIndex)*1e6;
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

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              MULTIPLE RUNS            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if false
    ets=16:20;
    Vos=[1e-1 1e0 1e1];
    outputM=zeros(length(ets)*length(Vos),3);
    tic
    
    %fault properties
    for j=1:length(Vos)
        evl.flt.Vo=0*evl.flt.l+Vos(j);
        
        for k=1:length(ets)
            % rheology
            evl.shz.n=0*evl.shz.n+3;
            evl.shz.m=0*evl.shz.n+2.7;
            evl.shz.etaM=evl.shz.etaM*0+10^ets(k)/(365*24*60*60)/1e6;
            evl.shz.etaK=evl.shz.etaM*0+10^ets(k)/(365*24*60*60)/1e6;
            
            s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-20e3)*sind(strike);
            d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)-20e3)*sind(strike);
            
            % pin patches with negative stress change at t=0
            flt.pinnedPosition=find((evl.evt{1}.KK{1,2}*(evl.evt{1}.src.slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evt{1}.rake)))<0);
            flt.pinnedPosition=sort(unique([flt.pinnedPosition; ...
                find(d<0)]));
            
            % initial conditions
            y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);
            % integration options
            options=ode.odeset('Refine',1,'RelTol',1e-6,'InitialStep',1e-9,'Stats','on','oDir','./postseismic_pinned');
            % provides evl.t and evl.y
            tpts=unique(sort([gpsdata{1}(:,1); gpsdata{2}(:,1); gpsdata{3}(:,1); gpsdata{4}(:,1)]));
            %tpts=[0:0.0027:0.7]
            evl.ode45(tpts,y0,options);
            
            gps.simulation(evl);
            datL2=zeros(length(gpsdata),1);
            for i=1:length(gpsdata)
                stationId=index{i};
                [C,ia,ib] = intersect(gpsdata{i}(end,1),tpts)
                modelN=gps.ua((stationId-1)*gps.vecsize+1,ib)+gps.uv((stationId-1)*gps.vecsize+1,ib);
                modelE=gps.ua((stationId-1)*gps.vecsize+2,ib)+gps.uv((stationId-1)*gps.vecsize+2,ib);
                modelD=gps.ua((stationId-1)*gps.vecsize+3,ib)+gps.uv((stationId-1)*gps.vecsize+3,ib);
                datL2(i)=sqrt((modelN-gpsdata{i}(end,2).*1e-3).^2+(modelE-gpsdata{i}(end,3).*1e-3).^2+(modelD-gpsdata{i}(end,4).*1e-3).^2)
            end
            outputM(k+length(ets)*(j-1),:)=[Vos(j) ets(k) sqrt(sum(datL2.^2))];
        end
    end
    toc
    p1=reshape(outputM(:,3),[5 3]);
    figure(99); hold on
    pcolor(log10(Vos),ets,p1)
    
end


%% GPS Generator

testxs=1e5*(-2.17:0.02:1.82)';
testys=1e5*(-6.38:0.02:2.22);
Nx=length(testxs);
Ny=length(testys);
testxs=repmat(testxs,1,Ny);
testys=repmat(testys,Nx,1);
polyx=japan_coasts_km(:,1)*1e3;
polyy=japan_coasts_km(:,2)*1e3;
inpoly=inpolygon(testxs,testys,polyx,polyy);
syngpsx=testxs(inpoly);
syngpsy=testys(inpoly);

return
fname='./gps/synthetic-gps.dat'
fid=fopen(fname,'wt');
fprintf(fid,'# export from unicycle\n');
fprintf(fid,'# n stationname x1 x2 x3\n');
fprintf(fid,'%05.5d %05.5d %f %f %f\n', ...
    [cumsum(ones(size(syngpsx,1),1)) cumsum(ones(size(syngpsx,1),1)) syngpsy syngpsx zeros(size(syngpsx))]');
fclose(fid);
%%
%syngps=[(1:size(syngpsx))' (1:size(syngpsx))' syngpsx syngpsy zeros(size(syngpsx))];
%save('./gps/synthetic-gps.dat','syngps','-ascii');
figure(5)
hold on
plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
plot(syngpsx,syngpsy,'b^','lineWidth',1,'markerfacecolor','b')
%plot(gps.x(:,1),gps.x(:,2),'r^','lineWidth',1)
%plot(gps_synthetic.x(:,2),gps_synthetic.x(:,1),'r^','lineWidth',1)
set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)

%%
gps_synthetic=unicycle.manifold.gpsReceiver('./gps/synthetic-gps.dat',evl,3);
%gps=unicycle.manifold.gpsReceiver('./gps/kumamoto-gps.dat',evl,3);

%% Plotting Simulated GPS
%gps_synthetic.simulation(evl);
  colormap default
    scale=2e7;
    figure(2001);clf;set(gcf,'Name','GPS synthetic maps')
    hold on; box on
   
    timeIndex=ceil(length(evl.t));
    %shz.plotShearZoneIndex();
    %shz.plotShearZoneEdges();
    %plot(gps.x(:,1),gps.x(:,2),'^');

     % plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
            axis equal

   %  scatter(gpsdata.y,gpsdata.x,200,gpsdata.dU,'filled')
     %     
   % scale*gps.ua(1:gps.vecsize:end,timeIndex)
    %evt{1}.plotPatch()
 %   toplot=gps_synthetic.ua(3:gps_synthetic.vecsize:end,timeIndex)+gps_synthetic.uv(3:gps_synthetic.vecsize:end,timeIndex);
    toplot=gps_synthetic.ua(3:gps_synthetic.vecsize:end,timeIndex);
    
  
    toplot=toplot+gps_synthetic.uv(3:gps_synthetic.vecsize:end,timeIndex);
          
    scatter(gps_synthetic.x(:,1),gps_synthetic.x(:,2),100,toplot,'filled','s')
    
   %  quiver(gpsdata.y,gpsdata.x, ...
    %     scale*gpsdata.dE,scale*gpsdata.dN,0,'r','lineWidth',1);
    % quiver(gps.x(:,1),gps.x(:,2), ...
     %       scale*(gps.ua(1:gps.vecsize:end,timeIndex)+gps.uv(1:gps.vecsize:end,timeIndex)), ...
      %      scale*(gps.ua(2:gps.vecsize:end,timeIndex)+gps.uv(2:gps.vecsize:end,timeIndex)),0,'b','lineWidth',2);
        

       
       % add in quiver command for data vectors as well as model
    
%     for k=1:length(gps.stationName)
%         text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
%     end
    
 %   stressComponent=3;
    %toplot=evl.evt{1}.KL{1,stressComponent}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,stressComponent}*(evt{1}.slip.*sind(evt{1}.rake));
    %shz.plotShearZoneEdges(toplot);
    %plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
            plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
                 plot(1e3*faults_maps2(:,1),1e3*faults_maps2(:,2),'r-','lineWidth',2)
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
   set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
    h=colorbar();
    ylabel(h,'Uplift (m)')
    title(sprintf('Surface displacements at time %d days',evl.t(timeIndex)/day));
    
    %%
   %clear all 
    load('./inversion_result.mat')
    

%% Uncertainty
       toplotnew=log10(1e6*J2(1:582))-log10((eI2./(86400*91)));
uncert=1e-13./(eI2./(86400*91));
weight=ones(size(uncert))*1e-5;
minmax(toplotnew)
prior=19;
visc=(weight*prior+1./uncert.*toplotnew)./(weight+1./uncert);
minmax(visc)
 figure(9) 
 clf
   hold on
   flt.plotPatch()
   %flt.plotPatch(rs)
   plot(1e3*volcano_km(:,1),1e3*volcano_km(:,2),'k^','markerfacecolor','y')
      shz{1}.plotShearZoneEdges(eI2(1:(end-6)))
  % shz{1}.plotShearZoneEdges(toplotnew(1:(end-6)))
  % shz{2}.plotShearZoneEdges(visc((end-5):end))
   %shz{3}.plotShearZoneEdges(ekk)
 %        shz{1}.plotShearZoneEdges(eI2(7:(end)))
  % shz{2}.plotShearZoneEdges(e11(1:6))
       colormap(jet);
            plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
                 plot(1e3*faults_maps2(:,1),1e3*faults_maps2(:,2),'r-','lineWidth',2)
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
     h=colorbar();
       set(gca,'clim',minmax((get(gca,'clim'))/1e0));
       set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)

%   shz{1}.plotShearZoneEdges(uncert(1:(end-6)))

    %%
   figure(9)
   hold on
   %flt.plotPatch(rs)
   toplotnew=log10(1e6*J2(1:582))-log10((eI2./(86400*91)));
      shz{1}.plotShearZoneEdges(uncert(1:(end-6)))
  % shz{1}.plotShearZoneEdges(toplotnew(1:(end-6)))
  % shz{2}.plotShearZoneEdges(toplotnew((end-5):end))
   %shz{3}.plotShearZoneEdges(ekk)
 %        shz{1}.plotShearZoneEdges(eI2(7:(end)))
  % shz{2}.plotShearZoneEdges(e11(1:6))
       colormap(jet);
            plot(1e3*japan_coasts_km(:,1),1e3*japan_coasts_km(:,2),'k-','lineWidth',1)
                 plot(1e3*faults_maps2(:,1),1e3*faults_maps2(:,2),'r-','lineWidth',2)
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    axis equal tight
     h=colorbar();
       set(gca,'xlim',[-180 180]*1e3,'ylim',[-200 140]*1e3,'zlim',[-100 0]*1e3)
       
       %%
       
       load '/Users/James/Dropbox/Kumamoto/Inversion/newresults/results_realdata_MCMC_w_outlier.mat'
       nz =[576 6 135];
       np = flt.N;

msig = sqrt(diag(Cm));

dhat = G*m;
gdhat = dhat(1:3*ngps);
ihat = dhat(3*ngps+1:end);

r = data - dhat;
grE = r(1:3:3*ngps);
grN = r(2:3:3*ngps);
grU = r(3:3:3*ngps);
ir = r(3*ngps+1:end);
gr = r(1:3*ngps);

rchi2 = r'/Cd*r ./ (length(data));
VR = (1 - (r'/Cd*r)/(data'/Cd*data));
iVR = 1 - sum((ir.^2./isig.^2))/sum((idata.^2./isig.^2));
gVR = 1 - sum((gr.^2./gsig.^2))/sum((gdata.^2./gsig.^2));

fprintf('# variance reduction. All: %10.6f%%, InSAR: %10.6f%%, GPS: %10.6f%%\n',VR*1e2,iVR*1e2,gVR*1e2);

ss = m(1:np);
ds = m(np+1:2*np);
rs = sqrt(ss.^2 + ds.^2);
e11 = m(2*np+0*(nz(1)+nz(2))+1:2*np+1*(nz(1)+nz(2)));
e12 = m(2*np+1*(nz(1)+nz(2))+1:2*np+2*(nz(1)+nz(2)));
e13 = m(2*np+2*(nz(1)+nz(2))+1:2*np+3*(nz(1)+nz(2)));
e22 = m(2*np+3*(nz(1)+nz(2))+1:2*np+4*(nz(1)+nz(2)));
e23 = m(2*np+4*(nz(1)+nz(2))+1:2*np+5*(nz(1)+nz(2)));
e33 = m(2*np+5*(nz(1)+nz(2))+1:2*np+6*(nz(1)+nz(2)));