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
% AUTHORS: 
% Sylvain Barbot and James Moore (July 29, 2016), Earth Observatory of Singapore

% clear all
% close all
% 
% if ispc
%     home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
% else
%     home = getenv('HOME');
% end
% eval(['addpath ' home '/Documents/src/unicycle/matlab'])
% 
% import unicycle.*
% 
% load './gmt/boundaries.dat';
% 
% gpsdata=cell(4,1);
% gpsdata{1} = importdata('./gps/KKN4out.dat');
% gpsdata{2} = importdata('./gps/NASTout.dat');
% gpsdata{3} = importdata('./gps/CHLMout.dat');
% gpsdata{4} = importdata('./gps/SNDLout.dat');
% 
% s2y=60*60*24*365;
% y2s=1./s2y;
% colors='kmcrgby';
% G=30e3;
% nu=1/4;
% minmax=@(x) [min(x),max(x)];
% heaviside=@(x) x>=0;
% omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);



function [dout]=postseismic_copula(minput)

import unicycle.*
load './gmt/boundaries.dat';
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

flt=geometry.triangleReceiver('./faults/qiu+15_1_receiver',greens.nikkhoo15(G,nu));
%flt=geometry.receiver('./faults/slab5',greens.okada92(G,nu));


% rate-and-state coefficients
flt.a=0*flt.l+minput(1);
% flt.b=flt.a-4e-3;

% reference velocity (m/yr)
%flt.Vo=0*flt.l+minput(2);
flt.Vo=0*flt.l+1e-6;

% characteristic slip distance
flt.l=0*flt.l+0.3;
flt.l=0*flt.l+7.5e-2;

% confining pressure (MPa)
flt.sigma=0*flt.l+1000;
flt.sigma=0*flt.l+1500;

% prevent afterslip in coseismic regions
flt.isRakeConstraint=true;
flt.rake=0*flt.l+90;

% OLD STUFF

% rate-and-state coefficients
flt.b=flt.a-8e-3;

strike=-72;
s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-0e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)+0e3)*sind(strike);
%index=find(1==omega(s/120e3).*omega(d/140e3));
index=find(1==omega(s/120e3).*omega(flt.xc(:,3)/40e3));
%flt.b(index)=flt.a(index)+4e-3;


% plate loading
flt.Vpl=0*flt.l+0.02;
azimuth=-160;
azimuthVector=[0*flt.l+sind(azimuth),0*flt.l+cosd(azimuth),0*flt.l];
flt.Vrake=atan2(sum(flt.dv.*azimuthVector,2),sum(flt.sv.*azimuthVector,2))*180/pi;
flt.rake=flt.Vrake;

% patches of interest
% flt.observationPoints={ ...
%     unicycle.geometry.observationPoint(0902), ... % trench
%     unicycle.geometry.observationPoint(2070), ... % Gorkha dÈcollement
%     unicycle.geometry.observationPoint(4688), ... % upper dÈcollement
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

evt={geometry.coseismicTriangle('./faults/qiu+15_1',0,greens.nikkhoo15(G,nu))};
%evt={geometry.coseismicPatch('./faults/test5.flt',0,greens.okada92(G,nu))};

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
%shz.n=0*shz.n+3;
%shz.m=0*shz.n+2.7;
shz.n=0*shz.n+1;
shz.m=0*shz.n+1;
shz.etaM=shz.etaM*0+minput(2)/(365*24*60*60)/1e6;
%shz.T0 = importdata('temps2.dat');
%shz.etaM=arr(shz.T0);

%shz.etaK=shz.etaM*0+1e16/(365*24*60*60)/1e6;
shz.etaK=1e-2*shz.etaM;

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
    %slip1=evl.evt{1}.src.slip;slip1(slip1<0.4)=NaN;
    %flt.plotPatch(slip);
    
    %sc=1e4;
    %flt.plotPatch(flt.b-flt.a);
    
    %flt.plotPatchIndex();
    %evt{1}.plotPatch(evt{1}.slip);
    
        shz.plotShearZoneEdges(log10(3.1536e+13*shz.etaK));
      %  shz.plotShearZoneEdges();
    %shz.plotUnitVectors(1e3);
    
    sc=1e4;
   % flt.plotPatch(flt.b-flt.a);
    flt.plotPatch();
    %flt.plotPatchIndex();
    evt{1}.plotPatch(evt{1}.slip);
  %  flt.plotPatch()
    
    %flt.plotPatch();
    
   % for k=1:length(flt.observationPoints)
    %    plot3( ...
    %        flt.observationPoints{k}.xc(1), ...
    %        flt.observationPoints{k}.xc(2), ...
    %        flt.observationPoints{k}.xc(3), ...
    %        [colors(k) '+'],'MarkerSize',100)
   % end
    
    plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    colormap jet
    h=colorbar('southoutside');
   
    ylabel(h,'Coseismic slip (m)');
    box on, grid on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
  %  title('(b-a) frictional properties')
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
    evl=ode.rateStrengtheningPower([],flt,shz,evt);
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

if false
    figure(1001);clf;set(gcf,'name','Stress interactions');
    subplot(2,1,1);gca;hold on
    
    %shz.plotShearZoneEdges();
    slip=evt{1}.slip;
    slip(7)=6;
    slip(1300)=6;
    slip(2436)=8;
    slip(85)=6;
    slip([2724,2776,1777,2509,1535,2274])=5.5;
    %toplot=evl.evt{1}.KK{1,2}*(slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(slip.*sind(evt{1}.rake));
    toplot=evl.evt{1}.KK{1,3}*(slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,3}*(slip.*sind(evt{1}.rake));
    %toplot(toplot>8)=8;
    %toplot=slip;
    
    %toplot=evl.KK{2,2};toplot=toplot(210,:);minmax(toplot)
    %flt.plotPatch(toplot);
    flt.plotPatch(toplot);
    %flt.plotPatchIndex();
    plot3(flt.xc(2436,1),flt.xc(2436,2),flt.xc(2436,3),'o','MarkerSize',10,'Color','b');
    title('Dip-slip traction component');
    sc=1e4;
    for k=1:length(flt.observationPoints)
        plot3( ...
            flt.observationPoints{k}.xc(1), ...
            flt.observationPoints{k}.xc(2), ...
            flt.observationPoints{k}.xc(3), ...
            [colors(k) '+'],'MarkerSize',100);
    end
    
    %plot(boundaries(:,1),boundaries(:,2),'k-','lineWidth',2)
    
    colormap(jet);
    h=colorbar();
    ylabel(h,'Stress');
    box on, grid on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    %set(gca,'view',[-227 42]);
    set(gca,'clim',[-1 1]*max((get(gca,'clim')))/1e0);
    axis equal
    %set(gca,'xlim',[-80 220]*1e3,'ylim',[-150 100]*1e3);
    fprintf('Stress interaction: review before simulation\n');
        subplot(2,1,2);cla;hold on;
    
    %toplot=evl.evt{1}.KL{1,5}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,5}*(evt{1}.slip.*sind(evt{1}.rake));
    %toplot=evl.evt{1}.KL{1,6}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,6}*(evt{1}.slip.*sind(evt{1}.rake));
    toplot=evl.evt{1}.KL{1,4}*(evt{1}.slip.*cosd(evt{1}.rake))+evl.evt{1}.KL{2,4}*(evt{1}.slip.*sind(evt{1}.rake));
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
%shz.etaM=arr(shz.T0);
%shz.etaK=1e-2*shz.etaM;
evl.shz.n=0*shz.n+3;
%evl.shz.n=0*shz.n+1;
evl.shz.m=0*shz.n+2.7;
%evl.shz.m=0*shz.n+1;
%evl.shz.etaM=shz.etaM*0+1e18/(365*24*60*60)/1e6;
%evl.shz.etaM=1e0*arr(evl.shz.T0);
%evl.shz.etaK=shz.etaM*0+1e17/(365*24*60*60)/1e6;
%evl.shz.etaK=1e0*arr(evl.shz.T0);
%evl.shz.etaK=1e-2*evl.shz.etaM;
evl.flt.Vo=0*evl.flt.l+3e0;
flt.l=0*flt.l+7.5e-2;
evl.flt.Vpl=0*flt.l+0.001;
evl.flt.b=evl.flt.a-8e-2;

s=+(flt.xc(:,2)+250e3)*cosd(strike)+(flt.xc(:,1)-20e3)*sind(strike);
d=+(flt.xc(:,1)+140e3)*cosd(strike)-(flt.xc(:,2)-20e3)*sind(strike);

% pin patches with negative stress change at t=0
flt.pinnedPosition=find((evl.evt{1}.KK{1,2}*(evl.evt{1}.src.slip.*cosd(evt{1}.rake))+evl.evt{1}.KK{2,2}*(evl.evt{1}.src.slip.*sind(evt{1}.rake)))<0);
flt.pinnedPosition=sort(unique([flt.pinnedPosition; ...
    find(d<0)]));

% initial conditions
y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);

% integration options
%options=ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','on','oDir','./postseismic_pinned');
options=ode.odeset('Refine',1,'RelTol',1e-3,'AbsTol',1e-4,'InitialStep',1e-6,'Stats','off');
% provides evl.t and evl.y

%tpts=[0:0.0027:0.7]
load 'tpts.mat';
evl.ode45(tpts,y0,options);

if ~exist('gps','var')
    gps=unicycle.manifold.gpsReceiver('./gps/gps_network.dat',evl,3);
%    gps.exportTimeSeriesToNED(evl.t','temp')
else
    gps.simulation(evl);
end


[~,stationName,~,~,~]=textread('./gps/gps_network.dat','%d %s %f %f %f','commentstyle','shell');

% ----------------------------------2------------------------------------------
% Load the observed data from the files and prep the output.
tst=[]; dout = [];
for i = 1:length(stationName)
    u1=0;
    u2=0;
    u3=0;
    filename=strcat('./gps/',strjoin(stationName(i)),'.ned');
    fid=fopen(filename,'r');
    out=textscan(fid,'%f %f %f %f %f %f %f','commentstyle','#');
    time=out{1};
    indexj=[];
    for j = 1:length(time)
        indexj=[indexj; find(tpts==time(j))];
    end
    
    u1=u1+gps.ua((i-1)*gps.vecsize+2,indexj)'+gps.uv((i-1)*gps.vecsize+2,indexj)';
    u2=u2+gps.ua((i-1)*gps.vecsize+1,indexj)'+gps.uv((i-1)*gps.vecsize+1,indexj)';
    u3=u3-gps.ua((i-1)*gps.vecsize+3,indexj)'-gps.uv((i-1)*gps.vecsize+3,indexj)';
    
    dout=[dout;u1];
    dout=[dout;u2];
    dout=[dout;u3];
    fclose(fid);
   
end


      
                
               
        
%dout = gps.exportTimeseriesArray();

end