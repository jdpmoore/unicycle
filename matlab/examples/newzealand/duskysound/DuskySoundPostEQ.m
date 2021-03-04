%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
import unicycle.*

%%%%%%%%%%%%% G: rigidity; nu: Possion's ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%
s2y = 60*60*24*365;y2s = 1./s2y;G = 30e3;nu = 1/4;
heaviside = @(x) x >= 0; 
omega  = @(x) heaviside(x+0.5)-heaviside(x-0.5);
minmax = @(x) [min(x),max(x)];

%%%%%%%%%%%% tectonic setting related data files %%%%%%%%%%%%%%%%%%%%%%%%%%
coast  = load('./auxdata/NZcoastlines_NEU.dat');
trench = load('./auxdata/DuskySound_trench.xy');


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

CoModel        = 'CP1_s23.5';
ReceiverFlt    = sprintf('./faults/beavan+10_%s_receiver0.flt',CoModel);
extendflt_d    = sprintf('./faults/beavan+10_%s_receiver_d.flt',CoModel);
extendflt_l    = sprintf('./faults/beavan+10_%s_receiver_l.flt',CoModel);
extendflt_r    = sprintf('./faults/beavan+10_%s_receiver_r.flt',CoModel);

earthModel_co  = greens.okada92(G,nu);
earthModel_shz = greens.shearZone16(G,nu);


flt = geometry.receiver({ReceiverFlt,extendflt_r,extendflt_l,extendflt_d},earthModel_co);

%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flt.a     = 0*flt.a     + 9e-4;     % rate-and-state coefficients
flt.Vo    = 0*flt.Vo    + 0.9;     % reference velocity (m/yr)
flt.sigma = 0*flt.sigma + 500;      % confining pressure (MPa)
%%%%%%%%%%% prevent afterslip in coseismic regions %%%%%%%%%%%%%%%%%%%%%%%%
flt.isRakeConstraint = true;
flt.Vrake            = 0*flt.Vrake + 90;



%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
CoSlip = sprintf('faults/beavan+10_%s.flt',CoModel);
evt    = {geometry.coseismicPatch(CoSlip,0,earthModel_co)};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S H E A R   Z O N E S                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
shzmodel = './faults/shearzone2';
shz      = geometry.shearZoneReceiver(shzmodel,earthModel_shz);

%%%%%%%%%%%%%%%% rheology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etaM0 = zeros(5,1);
etaM0(1:3) =  1e19;
etaM0(4:5) =  1e22;
shz.etaM   = CreateDepthDependent_viscovity(strcat(shzmodel,'_lvl.shz'),etaM0);
shz.etaM = shz.etaM/(s2y)/1e6;  
shz.etaK = shz.etaM;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           G E O M E T R Y   C H E C K                %
%                                                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


if false
    figure(1000);clf;set(gcf,'name','geometry');hold on
    shz.plotShearZoneEdges();
    flt.plotPatch();
    plot(coast(:,2),coast(:,1),'b')
    plot(trench(:,2),trench(:,1),'r.-')
    
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip.*sind(evt{k}.rake));
    end
    
    h=colorbar('southoutside');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic slip distribution')
    set(gca,'view',[47 40]);
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
     set(gca,'xlim',[-4.5e5,4.5e5],'ylim',[-4.5e5 4.5e5])
    return
end
%%%%%%%%%%%%%%%%%%%% viscosity setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
    figure(1001);clf;set(gcf,'name','shear zone and viscocity');hold on
    shz.plotShearZoneEdges(log10(shz.etaM*s2y*1e6));
    shz.plotShearZoneEdges();
    flt.plotPatch();
    plot(coast(:,2),coast(:,1),'k')
    plot(trench(:,2),trench(:,1),'r.-')
    colormap(jet)
    h=colorbar('southoutside');
    cmin = 16; cmax = 23;
    caxis([cmin cmax])
    
    ylabel(h,'viscosity');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic slip distribution')
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
    axis equal 
    set(gca,'xlim',[-4.5e5,4.5e5],'ylim',[-4.5e5 4.5e5])
    
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if ~exist('evl','var')
    tic
    % build stress kernels for integral equation
    evl = ode.rateStrengtheningBurgers([],flt,shz,evt,'./kernels/'); 
    toc
else
    % skip building stress kernel and update velocities & friction properties
    evl.flt = flt; 
    evl.shz = shz;
    evl.flt.dgf = 4; 
    evl.shz.dgf = 18;  
    evl.evt{1}.src = evt{1};
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S T R E S S   C H E C K                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false
    figure(1002);clf;hold on;
    % coseismic stress change
    ss = evt{1}.slip.*cosd(evt{1}.rake);
    ds = evt{1}.slip.*sind(evt{1}.rake);
    Ks = evl.evt{1}.KK{1,1}*ss+evl.evt{1}.KK{2,1}*ds;   %%% strike direction
    Kd = evl.evt{1}.KK{1,2}*ss+evl.evt{1}.KK{2,2}*ds;   %%% dip direction
    Kn = evl.evt{1}.KK{1,3}*ss+evl.evt{1}.KK{2,3}*ds;   %%% normal direction
    tau=(Ks.*cosd(flt.Vrake)+Kd.*sind(flt.Vrake));
    toplot = tau;
    
    flt.plotPatch(toplot);

    plot(coast(:,2),coast(:,1),'k')
    plot(trench(:,2),trench(:,1),'r.-')
   
    h=colorbar('southoutside');
    ylabel(h,'stress change');
    colormap(jet)
    box on;
        
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');

    axis equal 
    set(gca,'xlim',[-2.5e5,2.5e5],'ylim',[-2.5e5 2.5e5])
    
    figure(1003);clf;hold on
    flt.plotPatch();
    L11=evl.evt{1}.KL{1,1}*ss+evl.evt{1}.KL{2,1}*ds;
    L12=evl.evt{1}.KL{1,2}*ss+evl.evt{1}.KL{2,2}*ds;
    L13=evl.evt{1}.KL{1,3}*ss+evl.evt{1}.KL{2,3}*ds;
    L22=evl.evt{1}.KL{1,4}*ss+evl.evt{1}.KL{2,4}*ds;
    L23=evl.evt{1}.KL{1,5}*ss+evl.evt{1}.KL{2,5}*ds;
    L33=evl.evt{1}.KL{1,6}*ss+evl.evt{1}.KL{2,6}*ds;
    toplot = L13;
    
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    
    h=colorbar('southoutside');
    ylabel(h,'shear zone stress');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
 
    set(gca,'view',[47 40]);
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                 S I M U L A T I O N                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% pin patches with negative stress change at t=0
% coseismic stress change
 ss = evt{1}.slip.*cosd(evt{1}.rake);
 ds = evt{1}.slip.*sind(evt{1}.rake);
 Ks = evl.evt{1}.KK{1,1}*ss+evl.evt{1}.KK{2,1}*ds;   %%% strike direction
 Kd = evl.evt{1}.KK{1,2}*ss+evl.evt{1}.KK{2,2}*ds;   %%% dip direction
 Kn = evl.evt{1}.KK{1,3}*ss+evl.evt{1}.KK{2,3}*ds;   %%% normal direction
 tau=(Ks.*cosd(flt.Vrake)+Kd.*sind(flt.Vrake));
 
 pin1 = find(tau < 0); pin2 = find(evt{1}.slip > 0.7);
 pin  = unique([pin1;pin2]);
 
 flt.pinnedPosition = pin;

% initial condition
y0=zeros(1,flt.N*evl.flt.dgf+shz.N*evl.shz.dgf);
tic
% integration options
options=ode.odeset('Refine',1,'RelTol',1e-10,'InitialStep',1e-9,'Stats','on');
tsimu = 0:1/365.25:7;
evl.ode45(tsimu,y0,options);
toc

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
        gps=unicycle.manifold.gpsReceiver('./gps/gps_network.dat',evl,3);
    else
        gps.simulation(evl);
    end
%%%%%%%%%%%% export model data %%%%%%%%%%%%%%%%%%%%%%%%
system('mkdir -p ned_model');
gps.exportTimeSeriesToNED(evl.t','./ned_model')
end



%% GPS in map view

if true
    fid_dua = fopen('./ned_model/ua.dips.txt','w');
    fid_duv = fopen('./ned_model/uv.dips.txt','w');
    
    scale = 3e6;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    subplot(2,1,1);cla;hold on; box on;
    timeIndex=ceil(length(evl.t));
    plot(gps.x(:,1),gps.x(:,2),'^');
    toplot= sqrt(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2+evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2);
    flt.plotPatch(toplot);
    
%     shz.plotShearZoneEdges();
    plot(trench(:,2),trench(:,1),'r.-')
    h=colorbar();
    ylabel(h,'Slip (m)')
    
    quiver(0,200e3,scale*30e-3,0,0);
    for k=1:length(gps.stationName)
%         text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.ua((k-1)*gps.vecsize+1,timeIndex),scale*gps.ua((k-1)*gps.vecsize+2,timeIndex),0);
        fprintf(fid_dua,'%s  %15.8f  %15.8f   %10.6f    %10.6f  %10.6f\n',gps.stationName{k},...
            gps.x(k,1),gps.x(k,2),...
            gps.ua((k-1)*gps.vecsize+1,timeIndex),...
            gps.ua((k-1)*gps.vecsize+2,timeIndex),...
            gps.ua((k-1)*gps.vecsize+3,timeIndex));
    end
    plot(coast(:,2),coast(:,1),'k');
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
%     axis equal tight
    
    title(sprintf('Displacements from afterslip at time %f yr',evl.t(timeIndex)));
    axis equal 
    set(gca,'xlim',[-2e5,4e5],'ylim',[-2e5 3e5])
    
    subplot(2,1,2);cla;hold on; box on;
%     timeIndex=ceil(length(evl.t)/2);
    plot(gps.x(:,1),gps.x(:,2),'^');
%     flt.plotPatch();
     
    plot(trench(:,2),trench(:,1),'r.-')
    toplot=evl.y((evl.flt.N*evl.flt.dgf)+6:evl.shz.dgf:end,timeIndex);
    flt.plotPatch();
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
    plot(coast(:,2),coast(:,1),'k')
     
%     h=colorbar();
%     ylabel(h,'Strain')
    
    quiver(0,200e3,scale*30e-3,0,0);
    for k=1:length(gps.stationName)
        %text(gps.x(k,1)+2e3,gps.x(k,2)+2e3,gps.stationName{k})
        quiver(gps.x(k,1),gps.x(k,2),scale*gps.uv((k-1)*gps.vecsize+1,timeIndex),scale*gps.uv((k-1)*gps.vecsize+2,timeIndex),0);
        fprintf(fid_duv,'%s  %15.8f  %15.8f  %10.6f    %10.6f  %10.6f\n',gps.stationName{k},...
            gps.x(k,1),gps.x(k,2),...
            gps.uv((k-1)*gps.vecsize+1,timeIndex),...
            gps.uv((k-1)*gps.vecsize+2,timeIndex),...
            gps.uv((k-1)*gps.vecsize+3,timeIndex));
    end
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
%     axis equal tight
    title(sprintf('Displacements from viscoelastic flow at time %f yr',evl.t(timeIndex)));
    colorbar
    axis equal 
    set(gca,'xlim',[-2e5,4e5],'ylim',[-2e5 3e5])
    fclose(fid_dua);
    fclose(fid_duv);
end





%% 3d fault slip plot

if true
    
    
    figure(7);clf;
    set(gcf,'name','Source and receiver fault geometry')
    timeIndex=length(evl.t);
    
    subplot(2,1,1);gca;hold on;
    toplot=sqrt(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2+evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex).^2);
    flt.plotPatch(toplot);
%     flt.plotSlipVectors(evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),timeIndex),4e4)
    
    h=colorbar();
    ylabel(h,'slip (m)')
    
    title(sprintf('Afterslip (m) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;
    view([-30,20]);
    clim=get(gca,'clim');
%     set(gca,'clim',[0.0 0.045]);
    colormap(jet)
    box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
     axis equal 
%     set(gca,'xlim',[-2e5,4e5],'ylim',[-2e5 3e5])
    
    subplot(2,1,2);gca;hold on;
    
    toplot=evl.y((evl.flt.N*evl.flt.dgf)+3:evl.shz.dgf:end,timeIndex);
    shz.plotShearZoneEdges(toplot);
    shz.plotShearZoneEdges();
%     flt.plotPatch();
    
    h=colorbar();
    ylabel(h,'Strain')
    
    title(sprintf('Strain) t=%2.1e/%2.1e yr',evl.t(timeIndex),evl.t(end)));
    axis equal tight;
    view([-30,20]);
  
    clim=get(gca,'clim');
%     set(gca,'clim',[-1 1]*max(abs(clim(:))));
    box on, view([53 56]), xlabel('East (km)'), ylabel('North (km)')
     axis equal 
    set(gca,'xlim',[-2e5,4e5],'ylim',[-2e5 3e5])
end




