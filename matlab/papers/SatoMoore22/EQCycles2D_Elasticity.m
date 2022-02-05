% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         Earthquake cycles in 2D with both             %
%         on fault and off fault deformation            %
%                                                       %
% DESCRIPTION:                                          %
% Simulates unified earthquake cycles for slip on       %
% on a 2D dip-slip OR strike-slip fault                 %
%                                                       %
% OPTIONS:                                              %
% Off fault viscoelasticity                             %
% Surface gravity approximation                         %
% Time and space dependent parameters                   %
%                                                       %
% AUTHOR:                                               %
% James D. P. Moore (June, 2019)                        %
% Earth Observatory of Singapore                        %
% earth@jamesdpmoore.com                                %
%                                                       %
% REFERENCE:                                            %
% Sato D and Moore J D P                                %
% Displacements and stress associated with localised    %
% and distributed inelastic deformation with            %
% piecewise-constant elastic variations (GJI, 2022)     %
%                                                       %
% Please also run plotter.citations() at the bottom to  %
% generate appropriate additional citations.            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all
%close all
if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

minmax=@(x) [min(x(:)),max(x(:))]; % Minimum and maximum function
year = 60*60*24*365;
G = 30e3; % Shear modules in MPa
nu = 0.25; % Poisson's Ratio
g = 9.81; % Acceleration due to gravity

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                       OPTIONS                        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
tmax = 250*year; % Length of simulation
style = 1; % choose strike-slip (1) or dip-slip (2)
viscoelastic = false; % switch on and off viscoelastic regions
rheology = 1; % choose between maxwell (1) or power-law burgers (2)
gravity = false; % switch on and off gravitational effects
tstresses = false; % allow for time varying fault normal stresses
pinning = false; % allow for pinned fault patches
restart = false; % allow for continuing simulation from a previous one

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        FAULT                         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
flt = geometry.receiver({'./Faults/fault_90.seg'},greens.okada92(G,nu));
downdip = sqrt(flt.xc(:,1).^2+flt.xc(:,3).^2)./1e3;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    STRAIN REGIONS                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
if viscoelastic
    shz=geometry.shearZoneReceiver({'./faults/lowercrust.shz'},greens.shearZone16(G,nu));
    shz.x(:,1)=shz.x(:,1)+1;
    shz.xc(:,1)=shz.xc(:,1)+1;
    % rheology
    etaM=1e18; %Maxwell viscosity
    viscfac=1e6; %Scaling factor because the simulation tracks stresses in MPa
    shz.n=0*shz.n+1; %Maxwell element power exponent
    shz.m=0*shz.n+1; %Kelving element power exponent
    shz.etaM=shz.etaM*0+etaM/viscfac; 
    shz.etaK=shz.etaM./20; %Setting Kelvin viscosity equal to maxwell/20 (only relevent for Burgers)
else
    shz=[];
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   ELASTIC HETEROGENEITY              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
elast=unicycle.geometry.elasticEigendisplacement({'./faults/damage_zone.seg'}, unicycle.greens.okada92(G,nu));
% eKT = elast.tractionKernels(elast);
% eKD = elast.eigenDisplacementKernels(elast);
% [eKs,eKd]=unicycle.greens.computeDisplacementKernelsOkada85(elast,nu,elast.xc,3);
% xlocs=topo.x;
% xloc=topo.x(:,1);
% xlocs(:,2)=0*xlocs(:,2);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%        SURFACE GRAVITY AND OBSERVATION POINTS        %
%    Future version to include topographic effects     %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %
topo=geometry.shearZoneReceiver({'./faults/topography.shz'},greens.shearZone16(G,nu));
xlocs=topo.x;
xloc=topo.x(:,1);
xlocs(:,2)=0*xlocs(:,2);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%          velocity-weakening, a-b < 0                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
flt.sigma = flt.sigma*0 + 100; % effective confining pressure on fault (MPa)
%flt.sigma = flt.sigma*0 + 50 + 2*downdip;
flt.a = flt.a*0 + 1e-2;
flt.b = flt.a.*1.5;
flt.mu0 = flt.mu0*0 + 0.6; % static friction coefficient
flt.l = flt.l*0 + 0.012; % characteristic weakening distance (m)
flt.Vo = flt.Vo*0 + 1e-6; % reference slip rate (m/s)

% Plate loading
flt.Vs = flt.Vs*0 + 3e3; % shear wave speed (m/s)
flt.Vpl = flt.Vo*0 + 1e-9; % plate velocity (m/s)
if style == 2
    flt.Vpl = flt.Vpl.*cosd(flt.dip); % update plate velocity for dip-slip
end

% Velocity-strengthening at top and bottom ( a-b > 0 )
% 5-15km DOWN-DIP velocity-weakening
upper = downdip<5;
lower = downdip>15;
flt.b(upper) = flt.a(upper).*0.5;
flt.b(lower) = flt.a(lower).*0.5;

% Fault pinning
if pinning
    flt.pinnedPosition = find((downdip<15).*(downdip>5));   
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   Time varying stress                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if tstresses
    frequency = 1; %stress frequency wavelenth in years
    amplitude = 10; %-(downdip/4); %stress variation in MPa with spatial variation
    stressing=@(t) amplitude*sin(t./year./frequency);
else
    stressing=@(t) 0.*(downdip*t);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            G E O M E T R Y   C H E C K               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if true
    figure(1001);clf;set(gcf,'Color','White','name','geometry & properties');
    gca; hold on
    toplot = flt.a - flt.b;
    toplot(flt.pinnedPosition) = NaN;
    flt.plotPatch(toplot)
%     topo.plotShearZoneEdges();
    elast.plotPatch()
    if viscoelastic
        shz.plotShearZoneEdges(shz.etaM);
        shz.plotShearZoneEdges();
    end
    colormap jet
    h=colorbar('southoutside');
    ylabel(h,'Friction parameter (a-b)');
    box on, grid on;
    title('Friction Parameters');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[30 30]);
    if length(unique(minmax(toplot))) > 1
        set(gca,'clim',minmax(toplot));
    else
        set(gca,'clim',[0.9*min(toplot) 1.1*min(toplot)])
    end
    return
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%          S T R E S S   I N T E R A C T I O N         %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions if necessary
if ~exist('evl','var')
    % build stress kernels for integral equation
    tic
     evl=unicycle.ode.rateStrengtheningPower([],flt,shz,[], elast, './kernels/');
%     evl=unicycle.ode.rateStrengtheningPower([],flt,shz,[], './kernels/');
    if style == 1
        evl.shz.dgf=2+2*rheology;
    elseif style == 2
        evl.shz.dgf=3+3*rheology;
    end
    if ~viscoelastic
        evl.shz.N=0;
    end
    if gravity
        if ~exist(strcat(evl.knlpath,'GO_33.grd'),'file')
            G33=computeDisplacementKernelsSurfaceGravity(topo,nu,xlocs);
            grdwrite([0 1], [0 1], G33, [evl.knlpath,'GO_33.grd'])
        else
            [~,~,G33]=grdread([evl.knlpath,'GO_33.grd']);
        end
    else
        G33=[];
    end
    toc
else
    evl.flt=flt;
    evl.flt.dgf=4;
    if style == 1
        evl.shz.dgf=2+2*rheology;
    elseif style == 2
        evl.shz.dgf=3+3*rheology;
    end
    if ~viscoelastic
        evl.shz.N=0;
    end
end

if true
    figure(1002);clf;set(gcf,'Color','White','name','Stress interactions');
    gca; hold on
    slip=1-(upper+lower);
    toplot=evl.KK{style,style}*slip;
    toplot(flt.pinnedPosition) = NaN;
    flt.plotPatch(toplot)
    colormap jet
    h=colorbar('southoutside');
    ylabel(h,'Traction (MPa)');
    box on, grid on;
    if (style==1)
        title('Strike-slip traction component');
    else
        title('Dip-slip traction component');
    end
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[30 30]);
    %set(gca,'clim',minmax(toplot));
    fprintf('Stress interaction: review before simulation\n');
    return
end

%%

Kuu = evl.ED{1}(2:3:end,:);
KTu = evl.ET{1};
toInv = [Kuu -Kuu; KTu -KTu];
A = inv(toInv);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
if viscoelastic
    % choose between maxwell (1) or power-law burgers (2) - (Burgers NOT
    % currently working, will be fixed in the next version
    if rheology == 1
        yp=@(t,y) ode_acceleration_max(t,y,evl,style,stressing);
    elseif rheology == 2
        yp=@(t,y) ode_acceleration_pwr(t,y,evl,style,stressing);
    end
    y0=zeros(1,evl.flt.N*evl.flt.dgf+evl.shz.N*evl.shz.dgf);
else  
    yp=@(t,y) ode_acceleration_flt(t,y,evl,style,stressing);    
    y0=zeros(1,evl.flt.N*evl.flt.dgf);
end

% initial conditions
y0(3:evl.flt.dgf:(evl.flt.N*evl.flt.dgf)) = log(evl.flt.Vo./evl.flt.Vpl);
y0(4:evl.flt.dgf:(evl.flt.N*evl.flt.dgf)) = log(evl.flt.Vpl*0.98./evl.flt.Vo);
tstart = 0;

% True to resume simulation from previous saved conditions
storedstate = [evl.knlpath 'restart.mat'];
if restart
    if exist(storedstate)
        load(storedstate)
        if length(yend)==(evl.flt.N*evl.flt.dgf)
            y0 = yend;
            tstart = tend;
            tmax = tmax + tend;
        end
    end
end

% Solve the system
tic
options=odeset('Refine',1,'RelTol',1e-6,'AbsTol',1e-7,'InitialStep',1e-3,'MaxStep',year);
[evl.t,evl.y]=ode45(yp,[tstart tmax],y0,options);
toc
if restart
    % Save solution for restarting
    yend=evl.y(end,:);
    tend=evl.t(end);
    save(storedstate,'yend','tend')
end

%% Generate plotting object, then choose which plots to generate or things to export
plotter = Plotting(style,downdip,evl,xlocs,gravity,G33);
max(plotter.Vmax)
plotter.slipDeficit(1)
plotter.surfaceDisplacement([5 6],5)
plotter.slipContours(4,[0.1, year/12, 5*year],false)
plotter.velocityTime(2,stressing)
plotter.velocityTimestep(3,stressing)
if viscoelastic
    plotter.viscoelastic(10,250)
    plotter.viscoelasticTime(11,1)
end
tic; plotter.velocityEric(666,5000); plotter.slipContoursOver([1, year, 10*year]); toc
plotter.citations([viscoelastic gravity pinning tstresses])
%plotter.exportGMT('filename')