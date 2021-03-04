% This is an afterslip only model and is a part of the UNICYCLE program. 
% 
% AUTHOR: 
% Dongju Peng, Mar 2017.
% Sagar Masuti, Mar 2017.
% -------------------------------------------------------------------------

function [dout]=afterslip(minput)

G=30e3;
nu=1/4;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

flt=geometry.receiver('./faults/receiver.seg',greens.okada92(G,nu));

% rate-and-state coefficients
flt.a=0*flt.L+minput(1);
% flt.b=flt.a-4e-3;

% reference velocity (m/yr)
flt.Vo=0*flt.L+minput(2);

% characteristic slip distance
flt.l=0*flt.L+0.3;

% confining pressure (MPa)
flt.sigma=0*flt.L+1000;

% prevent afterslip in coseismic regions
flt.isRakeConstraint=true;
flt.Vrake=0*flt.L+90;

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

evt={geometry.coseismicPatch('faults/coseismic.flt',0,greens.okada92(G,nu))};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           G E O M E T R Y   C H E C K                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if false 
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    subplot(1,2,1);cla;hold on
    flt.plotPatch(flt.a-flt.b);
    
    h=colorbar('South');
    ylabel(h,'Friction properties');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('(a-b) frictional properties')
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(flt.b-flt.a)))
    axis equal
    
    subplot(1,2,2);cla;hold on
    flt.plotPatch();
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    
    h=colorbar('South');
    ylabel(h,'Slip (m)');
    box on;
    xlabel('East (m)');ylabel('North (m)'); zlabel('Up (m)');
    title('Coseismic slip distribution')
    set(gca,'view',[47 40]);
    axis equal
    fprintf('Sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl','var')
    % build stress kernels for integral equation
    evl=ode.ratestrengthening([],flt,[],evt,'./copula_kernels/');
else
    % skip building stress kernel and update velocities & friction properties
    evl.flt=flt;
    evl.shz = shz;
    evl.evt{1}.src=evt{1};
    evl.flt.dgf = 4; 
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                 S I M U L A T I O N                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% initial condition
y0=zeros(1,flt.N*evl.flt.dgf);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-10,'InitialStep',1e-9,'Verbose',false);
% creates variables evl.t and evl.y
tsimu=[0:0.01:2];
evl.ode45(tsimu,y0,options);


% GPS data (network.dat contains a list of stations with name and coordinates)
if ~exist('gps','var')
    gps=unicycle.manifold.gpsReceiver('./gps/gps_network.dat',evl,3);
%    gps.exportTimeSeriesToNED(evl.t','temp')
else
    gps.simulation(evl);
end

dout = gps.exportTimeseriesArray();


end

