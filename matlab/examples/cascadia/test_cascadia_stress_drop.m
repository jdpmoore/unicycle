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
%    [evl.t,evl.y]=evl.ode45((0:0.00274:2),y0,options);
%
% AUTHOR: 
% Sylvain Barbot (March 30, 2015), Earth Observatory of Singapore

clear all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

s2y=60*60*24*365;
y2s=1./s2y;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                E A R T H   M O D E L                 %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% triangle dislocations EGF
earthModelGimbutas12=greens.gimbutas12(30e3,1/4);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

rcv=geometry.triangleReceiver('./faults/schmidt_gao_Apr2001_receiver',earthModelGimbutas12);
%rcv=geometry.triangleReceiver('./faults/schmidt_gao_Aug2009_receiver',earthModelGimbutas12);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S O U R C E   F A U L T S               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

src=geometry.source();
%src=geometry.triangleReceiver('./faults/schmidt_gao_Apr2005_receiver',earthModelGimbutas12);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

evt={geometry.coseismicTriangle('./faults/schmidt+gao',0,earthModelGimbutas12)};
%evt={geometry.coseismicTriangle('./faults/schmidt_gao_Aug2009',0,earthModelGimbutas12)};


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S A N I T Y   C H E C K                %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

if true
    figure(1000);clf;set(gcf,'name','slip distribution & frictional properties');
    subplot(1,2,1);cla;hold on
    for k=1:length(evt)
        evt{k}.plotPatch(evt{k}.slip);
    end
    h=colorbar('South');
    title(h,'slip (m)');
    box on; grid on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    axis equal
    
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%             G R E E N S   F U N C T I O N S          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl')
    % build stress kernels for integral equation
    tic
    evl=ode.ratestrengthening(src,rcv,evt);
    toc
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.rcv=rcv;
    evl.evt=evt;
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%    C O S E I S M I C   S T R E S S   C H A N G E     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


if true
    tss=evl.evt{1}.Fss*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
        evl.evt{1}.Fds*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    tsd=evl.evt{1}.Fsd*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
        evl.evt{1}.Fdd*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake));
    
    stressDrop=-sum(tss.*(evl.evt{1}.src.slip.*cosd(evl.evt{1}.src.rake))+ ...
                    tsd.*(evl.evt{1}.src.slip.*sind(evl.evt{1}.src.rake)))/ ...
                sum(evl.evt{1}.src.slip);
    centroidDepth=-sum(evl.evt{1}.src.slip.*evl.evt{1}.src.xc(:,3))/ ...
                   sum(evl.evt{1}.src.slip);
           
    fprintf('stress drop: %f MPa\n',stressDrop);
    fprintf('centroid depth: %f km\n',centroidDepth/1e3);
    
    m=max(max(abs(tss)),max(abs(tsd)))/2;
    
    figure(2000);clf;set(gcf,'name','coseismic stress change');
    subplot(1,2,1);cla;hold on
    evt{1}.plotPatch(tss);
    
    h=colorbar('South');
    xlabel(h,'shear stress in the strike direction (MPa)');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*m);
    axis equal
    
    
    subplot(1,2,2);cla;hold on
    evt{1}.plotPatch(tsd);
    
    h=colorbar('South');
    title(h,'shear stress in the dip direction (MPa)');
    box on;
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    set(gca,'view',[47 40]);
    set(gca,'clim',[-1 1]*m);
    axis equal
    
    fprintf('coseismic stress change: review friction properties before simulation.\n');
    return
end


