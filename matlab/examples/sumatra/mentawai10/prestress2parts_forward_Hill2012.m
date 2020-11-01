% This script simulates postseismic deformation following the
% 2010 Mw 7.8 Mentawai tsunami earthquake with afterslip on
% the geometry that combines the planar geometry from Hill2012
% at shallower depth and Slab1.0 at deeper depth.
%
% Prestress conditions are specified separated for two regions.
% The southeastern region that has mostly slipped during the
% 2007 coseismic, aftershocks, and afterslip, while the 
% northwestern region that has no recent large earthquakes.
%
% Friction equilibrium is enforced on the "receiver" faults,
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
%    evt={geometry.coseismicPatch('./faults/coseismic.flt',0)};
%
% to obtain a afterslip distribution at given time steps (here daily 
% solutions for 2 years), use:
%
%    [evl.t,evl.y]=evl.ode45((0:0.00274:2),y0,options);
%
% AUTHORS: 
% Lujia Feng (May 2016), Earth Observatory of Singapore
%
% REFERENCE:
% Feng, L., S. Barbot, E. M. Hill, I. Hermawan, P. Banerjee, and D. H. Natawidjaja (2016), 
% Footprints of past earthquakes revealed in the afterslip of the 2010 Mw 7.8 Mentawai tsunami earthquake, 
% Geophys. Res. Lett., 43(18), 9518–9526, doi:10.1002/2016GL069870.

clear all
%close all

% change to your own path
addpath ~/SMT/SMT_MDL/EQ20101025/unicycle/matlab 
import unicycle.*

fltDir   = './faults/';
gpsDir   = './gps/';
plotDir  = './plots/';
figType  = '-depsc'; figExt = '.ps';   % for ps
figType  = '-dpng';  figExt = '.png'; % for png
figVisit = 'on';    % figure visibility

s2y=60*60*24*365;
y2s=1./s2y;

% convert from local to latlon
lon0   =  99.969549;
lat0   =  -4.229013;
strike = 322.406483;
rot    = strike-360;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% homogeneous elastic half space, Poisson's solid
% G=30e3;  shear modulus [MPa]
% nu=0.25; Poisson's solid
earthModel=greens.okada92(30e3,0.25);

time = 0:0.00274:5.001;

% read in coast lines
fcoastName = [gpsDir 'coasts.xyz'];
fin        = fopen(fcoastName,'r');
coastCell  = textscan(fin,'%f %f %f','CommentStyle','#');       % Delimiter default white space
coast      = cell2mat(coastCell);
coastyy    = coast(:,1);
coastxx    = coast(:,2);
coastzz    = coast(:,3);
ind = coastyy>4e5 | coastyy<-2e5;
coastyy(ind) = NaN;
coastxx(ind) = NaN;
coastzz(ind) = NaN;
fclose(fin);

% read in trench
ftrenchName = [gpsDir 'trench.xyz'];
fin        = fopen(ftrenchName,'r');
trenchCell  = textscan(fin,'%f %f %f','CommentStyle','#');      % Delimiter default white space
trench      = cell2mat(trenchCell);
trenchyy    = trench(:,1);
trenchxx    = trench(:,2);
trenchzz    = trench(:,3);
ind = trenchyy>=-2e5 & trenchyy<=4e5;
trenchyy = trenchyy(ind);
trenchxx = trenchxx(ind);
trenchzz = trenchzz(ind);
fclose(fin);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

rcv=geometry.receiver([fltDir 'MTW2010_receiver_rake90'],earthModel);

% depth
depth=-rcv.xc(:,3);
% rate-and-state coefficients
rcv.a=0*rcv.L+1e-2;
% characteristic slip distance
rcv.l=0*rcv.L+0.5; % does not matter for ratestrengthening

%-------------------------------------------------------
rcv.b=rcv.a-6e-4;       % misfit 2.6786
rcv.Vo=0*rcv.L+0.002;   % misfit 2.6786

theta2=2*(90-20+0*rcv.dip);

% depth
depth=-rcv.xc(:,3);
% constant confining pressure (MPa)
rcv.sigma=0*rcv.L+200; % rate-strengthening
rcv.tau=0*rcv.L+200;
rcv.mu0=0*rcv.L+1;
mind=480; % right 20 columns x 24 rows minus field
mind=456; % right 19 columns x 24 rows minus field
rind=rcv.id>mind;
%rind=rcv.id>480 & depth<=10e3; % right 18 columns x 24 rows minus field
rcv.tau(rind)=200.18;  % misfit 2.6786
lind=~rind;
rcv.tau(lind)=200.08;  % misfit 2.6786

% initial Coulomb stress change
cff0=(rcv.tau-rcv.mu0.*rcv.sigma);
rcv.Vpl=2*rcv.Vo.*sinh(cff0./((rcv.a-rcv.b).*rcv.sigma));
ind=cff0<0;
rcv.Vpl(ind)=0;

amb = rcv.a(1)-rcv.b(1);
Vo  = rcv.Vo(1);
lcff = rcv.tau(1)   - rcv.sigma(1);
rcff = rcv.tau(end) - rcv.sigma(1);

% rake constraint or not
rcv.isRakeConstraint=true;

options.sanityCheck=false;
if options.sanityCheck
    % (a-b)*sigma
    figh = figure('Visible',figVisit);
    clf; hold on;
    set(gcf,'name','geometry & frictional properties');
    toplot=(rcv.a-rcv.b).*rcv.sigma;
    rcv.plotPatch(toplot);
    plot3(coastxx,coastyy,coastzz,'k');
    h=colorbar; box on; ylabel(h,'(a-b)*sigma frictional properties');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    axis equal
    
    % depth
    figh = figure('Visible',figVisit);
    clf; hold on;
    set(gcf,'name','depth');
    rcv.plotPatch(rcv.xc(:,3)*1e-3);
    plot3(coastxx,coastyy,coastzz,'k');
    h=colorbar; box on; ylabel(h,'depth');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    axis equal
    
    % sigma
    figh = figure('Visible',figVisit);
    clf; hold on;
    set(gcf,'name','confining pressure');
    rcv.plotPatch(rcv.sigma);
    plot3(coastxx,coastyy,coastzz,'k');
    h=colorbar; box on; ylabel(h,'confining pressure');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    axis equal

    % cff0
    figh = figure('Visible',figVisit);
    clf; hold on;
    set(gcf,'name','pre-earthquake Coulomb stress');
    rcv.plotPatch(cff0);
    plot3(coastxx,coastyy,coastzz,'k');
    h=colorbar; box on; ylabel(h,'pre-earthquake Coulomb stress');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    axis equal
    return
end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%src=geometry.source([fltDir 'MTW2010_source_patch.flt']);
src=geometry.source();

%----------------------------------------- (1) base model -----------------------------------------

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
evt={};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% compute Green's functions only if necessary
evlMatfile = [fltDir 'evl_ratestrengthening_prestress.mat'];
% remove evl.mat file if geometry has changed!
%delete(evlMatfile);
% if evl exists in MATLAB Workspace
if exist('evl0')
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl0.Kss,1)==size(evl0.rcv.x,1),msg);
    assert(size(evl0.Fss,2)==size(evl0.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    % no need to update events if events have not changed!!!
    evl0.src=src;
    evl0.rcv=rcv;
    %evl.evt{1}.src.rake=evt{1}.rake;
% if evl Matfile exists
elseif exist(evlMatfile,'file')
    load(evlMatfile);
    evl0.src=src;
    evl0.rcv=rcv;
    %evl.evt{1}.src.rake=evt{1}.rake;
    fprintf('evl in %s is loaded.\n',evlMatfile);
else
    % build stress kernels for integral equation
    % full rate-and-state equations
    %evl=ode.rateandstate(src,rcv,evt);
    % rate-strengthening only
    evl0=ode.ratestrengthening_prestress(src,rcv,evt);
    save(evlMatfile,'evl0');
end

% integration options
%options=ode.odeset('Refine',1,'RelTol',1e-11,'AbsTol',1e-6,'InitialStep',1e-5);
options=ode.odeset('Refine',1,'RelTol',1e-8,'InitialStep',1e-5);
y0=zeros(1,rcv.N*evl0.dgf);
[sol0]=evl0.ode45(time,y0,options);

% GPS data (3 is vecsize)
gps0=unicycle.manifold.gpsReceiver([gpsDir 'EQ20101025po.sites'],evl0,3);
gps0.d=gps0.dataVectorFromTimeSeries(gpsDir,time); % filters data out of the range [time(1) time(end)]
d0=gps0.dataVectorFromModel(evl0);


%----------------------------------------- (2) afterslip model -----------------------------------------

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% ode45.m is the script where evts are added!!!
evt={geometry.coseismicPatch([fltDir 'MTW2010_coseismic_Hill2012.flt'],0,earthModel)};
%evt={};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% compute Green's functions only if necessary
evlMatfile = [fltDir 'evl_ratestrengthening_prestress_Hill2012.mat'];
% remove evl.mat file if geometry has changed!
%delete(evlMatfile);
% if evl exists in MATLAB Workspace
if exist('evl')
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    % no need to update events if events have not changed!!!
    evl.src=src;
    evl.rcv=rcv;
    %evl.evt{1}.src.rake=evt{1}.rake;
% if evl Matfile exists
elseif exist(evlMatfile,'file')
    load(evlMatfile);
    evl.src=src;
    evl.rcv=rcv;
    %evl.evt{1}.src.rake=evt{1}.rake;
    fprintf('evl in %s is loaded.\n',evlMatfile);
else
    % build stress kernels for integral equation
    % full rate-and-state equations
    %evl=ode.rateandstate(src,rcv,evt);
    % rate-strengthening only
    evl=ode.ratestrengthening_prestress(src,rcv,evt);
    save(evlMatfile,'evl');
end

options.sanityCheck=false;
if options.sanityCheck
    rak=evl.evt{1}.src.rake; % cosesimic rake
    slp=evl.evt{1}.src.slip; % coseismic slip

    % slip
    figh = figure('Visible',figVisit);
    clf;set(gcf,'name','coseismic event');
    rcv.plotPatch();
    %rcv.plotPatchIndex();
    for k=1:length(evt)
        fprintf('plotting slip distribution %s for time %f yr.\n',evt{k}.name,evt{k}.t0);
        evt{k}.plotPatch(evt{k}.slip);
    end
    h=colorbar; box on; ylabel(h,'coseismic slip (m)');
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    axis equal

    % coseismic stress change in dip direction
    figh = figure('Visible',figVisit);
    clf;set(gcf,'name','coseismic stress change');
    stress_dip=evl.evt{1}.Fsd*(slp.*cosd(rak))+...  % slip in strike direction
               evl.evt{1}.Fdd*(slp.*sind(rak));     % slip in dip direction
    rcv.plotPatch(stress_dip);
    h=colorbar(); box on;
    ylabel(h,'stress change in dip direction (MPa)')
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[60 40]);
    axis equal

    % coseismic stress change in strike direction
    figh = figure('Visible',figVisit);
    clf;set(gcf,'name','coseismic stress change');
    stress_str=evl.evt{1}.Fss*(slp.*cosd(rak))+...  % slip in strike direction
               evl.evt{1}.Fds*(slp.*sind(rak));     % slip in dip direction
    rcv.plotPatch(stress_str);
    h=colorbar(); box on;
    ylabel(h,'stress change in strike direction (MPa)')
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[60 40]);
    axis equal
    return
end

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-12,'AbsTol',1e-9,'InitialStep',1e-8,'Stats','on');
y0=zeros(1,rcv.N*evl.dgf);
[sol]=evl.ode45(time,y0,options);

% GPS data (network.dat contains a list of stations with name and coordinates)
gps=unicycle.manifold.gpsReceiver([gpsDir 'EQ20101025po.sites'],evl,3);
gps.d=gps.dataVectorFromTimeSeries(gpsDir,time); % filters data out of the range [time(1) time(end)]
dmodel=gps.dataVectorFromModel(evl);
dreal=dmodel-d0;

% misfit
misfit=sqrt(sum((gps.d-dreal).^2));

% subtract the base model from the afterslip model
gps.ur = gps.ur - gps0.ur;

% save time series predictions for stations
[ sites,locs ] = GPS_readsites('gps/SuGAr.sites');
fpo = fopen(['gps/EQ201010254029_Mw7.80po_Hill2012_prestress2parts_fwd_5yrs_mind' num2str(mind,'%3d') '_amb' num2str(amb,'%10.2e') '_Vo' num2str(Vo,'%10.2e') '_lcff' num2str(lcff,'%10.2e') '_rcff' num2str(rcff,'%10.2e') '.vel' ],'w');
fprintf(fpo,'# 1   2   3   4      5  6  7  8   9   10  11      12 13  14 15    16\n');
fprintf(fpo,'# Sta Lon Lat Height VE VN VU ErE ErN ErU Weight  T0 T0D T1 T1D   Cne\n');
fprintf(fpo,'# Height [m] Disp [mm] error [mm]\n');
for ii=1:length(gps.stationName)
   % one way
   statName = gps.stationName{ii};
   Mtt = evl.t';
   Mee = gps.ur((ii-1)*gps.vecsize+1,:)';
   Mnn = gps.ur((ii-1)*gps.vecsize+2,:)';
   Muu = gps.ur((ii-1)*gps.vecsize+3,:)';
   fsite = fopen(['gps/' statName '_Hill2012_prestress2parts_fwd_mind' num2str(mind,'%3d') '_amb' num2str(amb,'%10.2e') '_Vo' num2str(Vo,'%10.2e') '_lcff' num2str(lcff,'%10.2e') '_rcff' num2str(rcff,'%10.2e') '.out' ],'w');
   fprintf(fsite,'#time(yr) ENU (m) ENU_err (m)\n');
   enu = [ Mtt Mnn Mee -Muu ];
   fprintf(fsite,'%14.6f %12.5f %12.5f %12.5f\n',enu'); 
   fclose(fsite);

   % save po.vel
   ind = strcmpi(statName,sites);
   fprintf(fpo,'%4s %14.9f %14.9f %10.4f %7.1f %7.1f %7.1f %6.1f %4.1f %4.1f %6.1f\n',...
	   statName,locs(ind,:),Mee(end)*1e3,Mnn(end)*1e3,Muu(end)*1e3,0.0,0.0,0.0,1.0);
end
fclose(fpo);

% save accumulative slip after 5 years
fslip  = fopen(['results/Hill2012_prestress2parts_fwd_5yrs_mind' num2str(mind,'%3d') '_amb' num2str(amb,'%10.2e') '_Vo' num2str(Vo,'%10.2e') '_lcff' num2str(lcff,'%10.2e') '_rcff' num2str(rcff,'%10.2e') '.slip'],'w');
fslip0 = fopen(['results/Hill2012_prestress2parts_fwd_5yrs_mind' num2str(mind,'%3d') '_amb' num2str(amb,'%10.2e') '_Vo' num2str(Vo,'%10.2e') '_lcff' num2str(lcff,'%10.2e') '_rcff' num2str(rcff,'%10.2e') '_smallregion.slip'],'w');
%timeInd=find(abs(evl.t-3)<0.00137);
xx = rcv.xc(:,1); % center east
yy = rcv.xc(:,2); % center north
uu = rcv.xc(:,3); % center up
[ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
ss   = evl.y(1:evl.dgf:end,end);
ds   = evl.y(2:evl.dgf:end,end)-rcv.Vpl*evl.t(end);
ts   = zeros(size(ss));
slip = sqrt(ss.^2+ds.^2);
rake = atan2d(ds,ss);
[ ee,nn ] = rotate_xy(ss,ds,strike-360-90);
slipout = [ lon lat -uu*1e-3 ss ds ts slip rake ee nn ];
moment = sum(slip.*rcv.L.*rcv.W.*earthModel.G*1e6);
Mw     = log10(moment)*2/3-6.06667;
fprintf(fslip,'# (1)clon (2)clat (3)cdepth[km] (4)ss[m] (5)ds[m] (6)ts[m] (7)slip[m] (8)rake (9)east[m] (10)north[m]\n');
fprintf(fslip,'# moment = %10.5e Mw = %6.2f minslip = %10.4f [m] maxslip = %10.4f [m]\n',moment,Mw,min(slip),max(slip));
fprintf(fslip,'%14.6f %14.6f %8.2f %12.5f %12.5f %12.5f %12.5f %10.5f %12.5f %12.5f\n',slipout'); 
fclose(fslip);

kind = rcv.id>9*24 & rcv.id<=28*24 & rcv.x(:,3)>-20e3; % 24 rows minus field
slipout = slipout(kind,:);
fprintf(fslip0,'# (1)clon (2)clat (3)cdepth[km] (4)ss[m] (5)ds[m] (6)ts[m] (7)slip[m] (8)rake (9)east[m] (10)north[m]\n');
fprintf(fslip0,'# moment = %10.5e Mw = %6.2f minslip = %10.4f [m] maxslip = %10.4f [m]\n',moment,Mw,min(slip(kind)),max(slip(kind)));
fprintf(fslip0,'%14.6f %14.6f %8.2f %12.5f %12.5f %12.5f %12.5f %10.5f %12.5f %12.5f\n',slipout'); 
fclose(fslip0);

% interpolate data time series to model results 
pos=0; % position for d vector
gpsdur = zeros(size(gps.ur));
for ii=1:length(gps.stationName)
   % read data (modeled version) directly
   statName = gps.stationName{ii};
   finName  = [ 'gps/' statName '_clean_sunda_EQ20101025po_2model.rneu' ];
   [ rneu ] = GPS_readrneu(finName,1,1); % this matlab function is in ../mentawai10/gps/ folder
   Dtt  = rneu(:,2) - 2010.815068; % reference to EQ day
   Dnn  = rneu(:,3);
   Dee  = rneu(:,4);
   Duu  = rneu(:,5);

   iDnn = interp1(Dtt,Dnn,evl.t);
   iDee = interp1(Dtt,Dee,evl.t);
   iDuu = interp1(Dtt,Duu,evl.t);

   gpsdur((ii-1)*gps.vecsize+1,:) = iDee;
   gpsdur((ii-1)*gps.vecsize+2,:) = iDnn;
   gpsdur((ii-1)*gps.vecsize+3,:) = iDuu;
end

options.saveMovie = false
if options.saveMovie
    % save incremental slip & displacement for 5 years
    [ sites,locs ] = GPS_readsites('gps/SuGAr.sites');
    dirname = 'results/movie/';
    evlnum = length(evl.t);
    for kk=1:evlnum
         % save slip
         ss   = evl.y(1:evl.dgf:end,kk);
         ds   = evl.y(2:evl.dgf:end,kk)-rcv.Vpl*evl.t(kk);
         ts   = zeros(size(ss));
         slip = sqrt(ss.^2+ds.^2);
         rake = atan2d(ds,ss);
         [ ee,nn ] = rotate_xy(ss,ds,strike-360-90);
         slipout = [ lon lat -uu*1e-3 ss ds ts slip rake ee nn ];
         moment = sum(slip.*rcv.L.*rcv.W.*earthModel.G*1e6);
         Mw     = log10(moment)*2/3-6.06667;
         
         % save slip in smallregion
         kind = rcv.id>9*24 & rcv.id<=28*24 & rcv.x(:,3)>-20e3; % 24 rows minus field
         slipout = slipout(kind,:);
         foutName = [ dirname 'slip_evolution_' num2str(kk,'%06d') '.slip' ];
         fout = fopen(foutName,'w');
         fprintf(fout,'# (1)clon (2)clat (3)cdepth[km] (4)ss[m] (5)ds[m] (6)ts[m] (7)slip[m] (8)rake (9)east[m] (10)north[m]\n');
         fprintf(fout,'# moment = %10.5e Mw = %6.2f minslip = %10.4f [m] maxslip = %10.4f [m]\n',moment,Mw,min(slip),max(slip));
         fprintf(fout,'# time = %12.6f [yr]\n',evl.t(kk));
         fprintf(fout,'%14.6f %14.6f %8.2f %12.5f %12.5f %12.5f %12.5f %10.5f %12.5f %12.5f\n',slipout'); 
         fclose(fout);
    
         % save displacement prediction
         foutName = [ dirname 'displacement_predicted_' num2str(kk,'%06d') '.vel' ];
         fout = fopen(foutName,'w');
         fprintf(fout,'# 1   2   3   4      5  6  7  8   9   10  11      12 13  14 15    16\n');
         fprintf(fout,'# Sta Lon Lat Height VE VN VU ErE ErN ErU Weight  T0 T0D T1 T1D   Cne\n');
         fprintf(fout,'# Height [m] Disp [mm] error [mm]\n');
         fprintf(fout,'# time = %12.6f [yr]\n',evl.t(kk));
         for ii=1:length(gps.stationName)
            statName = gps.stationName{ii};
            Mee = gps.ur((ii-1)*gps.vecsize+1,kk)';
            Mnn = gps.ur((ii-1)*gps.vecsize+2,kk)';
            Muu = gps.ur((ii-1)*gps.vecsize+3,kk)';
            % save po.vel
            ind = strcmpi(statName,sites);
            fprintf(fout,'%4s %14.9f %14.9f %10.4f %7.1f %7.1f %7.1f %6.1f %4.1f %4.1f %6.1f\n',...
         	        statName,locs(ind,:),Mee*1e3,Mnn*1e3,Muu*1e3,0.0,0.0,0.0,1.0);
         end
         fclose(fout);
    
         % save displacement observation
         foutName = [ dirname 'displacement_observed_' num2str(kk,'%06d') '.vel' ];
         fout = fopen(foutName,'w');
         fprintf(fout,'# 1   2   3   4      5  6  7  8   9   10  11      12 13  14 15    16\n');
         fprintf(fout,'# Sta Lon Lat Height VE VN VU ErE ErN ErU Weight  T0 T0D T1 T1D   Cne\n');
         fprintf(fout,'# Height [m] Disp [mm] error [mm]\n');
         fprintf(fout,'# time = %12.6f [yr]\n',evl.t(kk));
         for ii=1:length(gps.stationName)
            statName = gps.stationName{ii};
            Dee = gpsdur((ii-1)*gps.vecsize+1,kk)';
            Dnn = gpsdur((ii-1)*gps.vecsize+2,kk)';
            Duu = gpsdur((ii-1)*gps.vecsize+3,kk)';
            % save po.vel
            ind = strcmpi(statName,sites);
            fprintf(fout,'%4s %14.9f %14.9f %10.4f %7.1f %7.1f %7.1f %6.1f %4.1f %4.1f %6.1f\n',...
         	        statName,locs(ind,:),Dee*1e3,Dnn*1e3,Duu*1e3,0.0,0.0,0.0,1.0);
         end
         fclose(fout);
    end
end

%% surface data simulation
options.plotGpsTimeSeries=true;
if options.plotGpsTimeSeries
    pos=0; % position for d vector
    for ii=1:length(gps.stationName)
       Dtt=gps.epochs{ii};
       Dnum=length(Dtt);
       Dee=gps.d(pos+1:pos+Dnum);
       Dnn=gps.d(pos+Dnum+1:pos+Dnum*2);
       Duu=gps.d(pos+Dnum*2+1:pos+Dnum*3);
       pos=pos+Dnum*gps.vecsize;

       figure;clf;set(gcf,'Name','GPS time series')
       statName=gps.stationName{ii};
       subplot(3,1,1);cla;
       hold on, box on;
       comp=1;
       plot(Dtt,Dee,'go');
       gps.plotTimeSeries(evl.t,ii,comp);
       title(sprintf('station %s, east displacement',statName));
       ylabel('east displacement (m)');
       
       subplot(3,1,2);cla;
       hold on, box on;
       comp=2;
       plot(Dtt,Dnn,'go');
       gps.plotTimeSeries(evl.t,ii,comp);
       title(sprintf('station %s, north displacement',statName));
       ylabel('north displacement (m)');
       
       subplot(3,1,3);cla;
       hold on, box on;
       comp=3;
       plot(Dtt,Duu,'go');
       gps.plotTimeSeries(evl.t,ii,comp);
       title(sprintf('station %s, vertical displacement',statName));
       ylabel('vertical displacement (m)'), xlabel('time (yr)');
    end
end

%% slip evolution movie
options.plotEvolution=true;
if options.plotEvolution

   figh = figure(1001); clf;
   ii=0;
   krange = [1:30:length(evl.t) length(evl.t)];
   knum = length(krange);
   for nn=1:knum
        kk   = krange(nn);
        figure(1001); cla; hold on; box on; caxis([0 2.2]);
        ss   = evl.y(1:evl.dgf:end,kk);
        ds   = evl.y(2:evl.dgf:end,kk)-rcv.Vpl*evl.t(kk);
        slip = sqrt(ss.^2+ds.^2);
        rake = atan2d(ds,ss);

        rcv.plotPatch(slip);
        plot3(coastxx,coastyy,coastzz,'k');
        plot3(trenchxx,trenchyy,trenchzz,'LineStyle','-','LineWidth',1,'Color','k','Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4);
        h = colorbar;
        title(sprintf('Slip [m] @ %2.1f years',evl.t(kk)));
        %axis equal tight; colorbar;
        set(gca,'xlim',[-3e5 3e5],'ylim',[-1e5 4e5],'zlim',[-50 0]*1e3);
        %clim=get(gca,'clim');
        %set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
        set(gca,'view',[60 30]);
        daspect([3 3 1]);
        figName = [ plotDir 'movie/slip_evolution_' num2str(ii,'%3d') figExt ];
        print(figh,figType,figName);
        
	ii=ii+1;
        pause(0.1)
   end

   component=2; % dip slip
   figure(1002); clf;
   for k=2:50:length(evl.t)
       figure(1002);
       subplot(2,1,1);cla;hold on
       toplot=evl.y(component:evl.dgf:end,k)-rcv.Vpl*evl.t(k);
       rcv.plotPatch(toplot);
       plot3(coastxx,coastyy,coastzz,'k');
       colorbar
       title(sprintf('dip slip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
       axis equal tight; colorbar;
       set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-165 0]*1e3);
       clim=get(gca,'clim');
       set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
       box on, view([41 30]), xlabel('x (km)'), ylabel('y (km)')
       
       subplot(2,1,2);cla;hold on
       toplot=((evl.y(component:evl.dgf:end,k)-evl.y(component:evl.dgf:end,k-1))/(evl.t(k)-evl.t(k-1))-rcv.Vpl)*y2s;
       rcv.plotPatch(toplot);
       plot3(coastxx,coastyy,coastzz,'k');
       colorbar
       title(sprintf('dip slip velocity (m/s) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
       axis equal tight; colorbar;
       set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-165 0]*1e3);
       clim=get(gca,'clim');
       set(gca,'clim',[-12 -1]);
       box on, view([41 30]), xlabel('x (km)'), ylabel('y (km)')
       
       pause(0.0125)
   end
end
