% DESCRIPTION:
% UNICYCLE script for the 2014 Mw 6.5 Napa Valley earthquake
% postseismic relaxation
%
% AUTHOR:
% Sylvain Barbot (Nov 24, 2014), Earth Observatory of Singapore

%clear all

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

% homogeneous elastic properties
G=30e3;
nu=1/4; % Poisson's solid

s2y=60*60*24*365;
y2s=1./s2y;

load 'share/coasts_km.dat'
load 'share/coasts_dim.dat'
load 'share/ca_faults_km.dat'
load 'share/ca_faults_dim.dat'

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

heaviside=@(x) x>=0;
omega=@(x) heaviside(x+0.5)-heaviside(x-0.5);

rcv=geometry.receiver({'./faults/receiver'});

% is rake constrained?
rcv.isRakeConstraint=true;

% % % % % % % % % % % % % % % % %
%
%  rate-and-state coefficients
%
% % % % % % % % % % % % % % % % %

% setup default properties for the Wharton Basin earthquake faults
try 
    pos=rcv.segments{1}.starti+(1:rcv.segments{1}.nPatch);
    % weakening properties
    rcv.a(pos)=1e-2;
    rcv.b(pos)=rcv.a(pos)-10e-3;
    % reference velocity (m/yr)
    rcv.Vo(pos)=10;
    % characteristic slip distance
    rcv.l(pos)=0.3;
    % confining pressure (MPa)
    rcv.sigma(pos)=600;
catch me
    if strcmp(me.identifier,'MATLAB:Containers:Map:NoKey')
        fprintf('could not find base name ''%s''; ignoring.\n',basename);
    else
        rethrow(me);
    end
end

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

evt={geometry.coseismicpatch('./faults/swei+14.flt',0)};

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           C O S E I S M I C   E V E N T S            %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

options.sanityCheck=false;
if options.sanityCheck
    xlim=[-20 20]*1e3;
    ylim=[-20 20]*1e3;
    figure(1000);clf;set(gcf,'name','geometry & frictional properties');
    hold on
    %rcv.plotPatch(rcv.Vo);
    %rcv.plotPatchIndex();
    for k=1:length(evt)
        fprintf('plotting slip distribution %s for time %f yr.\n',evt{k}.name,evt{k}.t0);
        evt{k}.plotPatch(evt{k}.slip);
    end
    down=2;
    scale=3e4;
    %quiver3(rcv.xc(1:down:end,1),rcv.xc(1:down:end,2),rcv.xc(1:down:end,3), ...
    %    scale*rcv.nv(1:down:end,1),scale*rcv.nv(1:down:end,2),scale*rcv.nv(1:down:end,3),0,'color','r');
    %quiver3(rcv.xc(1:down:end,1),rcv.xc(1:down:end,2),rcv.xc(1:down:end,3), ...
    %    scale*rcv.sv(1:down:end,1),scale*rcv.sv(1:down:end,2),scale*rcv.sv(1:down:end,3),0,'color','g');
    %quiver3(rcv.xc(1:down:end,1),rcv.xc(1:down:end,2),rcv.xc(1:down:end,3), ...
    %    scale*rcv.dv(1:down:end,1),scale*rcv.dv(1:down:end,2),scale*rcv.dv(1:down:end,3),0,'color','b');
    unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
    unicycle.plot.plot_faults(ca_faults_km*1e3,ca_faults_dim,xlim,ylim,[1 0.1 0.1]);
    h=colorbar(); box on
    ylabel(h,'Vo coseismic slip (m)')
    
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('coseismic slip')
    axis equal tight
    %set(gca,'view',[45 48],'clim',[0 1]*50,'xlim',xlim,'ylim',ylim);
    set(gca,'view',[-70 40],'clim',[0 1]*.6,'xlim',xlim,'ylim',ylim);
    fprintf('sanity check: review geometry and properties before simulation\n');
    return
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%              S L I P   E V O L U T I O N             %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% compute Green's functions only if necessary
if ~exist('evl')
    fprintf('build stress interaction kernels.\n');
    % build stress kernels for integral equation
    evl=ode.ratestrengthening(src,rcv,evt,G,nu);
    
    %
    % or try using the full rate-and-state equations
    %
    % evl=ode.rateandstate(src,rcv,evt,G,nu);
    %
    % and compare.
    
else
    msg='Source and/or receiver''s geometry have changed. clear all and re-run.';
    assert(size(evl.Kss,1)==size(evl.rcv.x,1),msg);
    assert(size(evl.Fss,2)==size(evl.src.x,1),msg);
    % skip building stress kernel and update velocities & friction properties
    evl.src=src;
    evl.rcv=rcv;
    %evl.cap=2;
    evl.evt{1}.src.rake=evt{1}.rake;
    %evl.evt=evl.eventKernels(evt,rcv,G,nu);
end

% GPS data (network.dat contains a list of stations with name and coordinates)
gps=unicycle.manifold.okada85.gps('./gps/opts.dat',evl,3);
%gps.d=gps.dataVectorFromTimeSeries('../gps/ps',evl.t);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%    C O S E I S M I C   S T R E S S   C H A N G E     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

options.stressChange=false;
if options.stressChange
    xlim=[-25 25]*1e3;
    ylim=[-25 25]*1e3;
    scale=2e4;
    figure(2000);clf;set(gcf,'name','coseismic stress change');
    rcv.plotPatch();
    %rcv.plotPatchIndex();
    for k=1:length(evt)
        fprintf('plotting coseismic stress change for event %s for time %f yr.\n',evt{k}.name,evt{k}.t0);
        %rcv.plotPatch(evl.evt{k}.Fsd*(evt{k}.slip.*cosd(evt{k}.rake))+...
        %              evl.evt{k}.Fdd*(evt{k}.slip.*sind(evt{k}.rake)));
        rcv.plotPatch(-(evl.evt{k}.Fss*(evt{k}.slip.*cosd(evt{k}.rake))+...
                        evl.evt{k}.Fds*(evt{k}.slip.*sind(evt{k}.rake))));
        evt{k}.plotPatch();
    end
    unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
    unicycle.plot.plot_faults(ca_faults_km*1e3,ca_faults_dim,xlim,ylim,[0.8 0.8 0.1]);
    
    t=repmat(evl.evt{k}.Fss*(evt{k}.slip.*cosd(evt{k}.rake))+...
             evl.evt{k}.Fds*(evt{k}.slip.*sind(evt{k}.rake)),1,3).*rcv.sv+...
      repmat(evl.evt{k}.Fsd*(evt{k}.slip.*cosd(evt{k}.rake))+...
             evl.evt{k}.Fdd*(evt{k}.slip.*sind(evt{k}.rake)),1,3).*rcv.dv;
    sscale=1e2;
    quiver3(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3),...
        sscale*t(:,1),sscale*t(:,2),sscale*t(:,3),0,'color',[0.1,0.3,0.3],'linewidth',2)
    
    for i=1:length(gps.stationName)
        [gt,gn,ge,gd,gsn,gse,gsd]=textread(['./gps/' gps.stationName{i} '.dat'],...
                '%f %f %f %f %f %f %f','commentstyle','shell');
        [~,index]=min(abs(gt-gt(1)-1.0));
        quiver(gps.x(i,1),gps.x(i,2),scale*ge(index),scale*gn(index),'k','linewidth',2);
        h=plot.ellipse(scale*0.01,scale*0.01,0,gps.x(i,1)+scale*ge(index),gps.x(i,2)+scale*gn(index),'k');
        set(h,'linewidth',2)
    end
    h=colorbar(); box on
    ylabel(h,'coseismic stress change (MPa)')
    
    xlabel('east (m)');ylabel('north (m)'); zlabel('up (m)');
    title('coseismic slip')
    axis equal
    set(gca,'view',[-70 40],'clim',[-1 1]*5,'xlim',xlim,'ylim',ylim);
    fprintf('review stress change before simulation\n');
    return
end

%% solver

time=(0:0.0027:0.100);

% initial condition
y0=zeros(1,rcv.N*evl.dgf);

% integration options
options=ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5);
[evl.t,evl.y]=evl.ode45(time,y0,options);

assert(0==sum(sum(isnan(evl.y))),'unicycle:invalid simulation')

%% export to Paraview

options.exportVtp=false;
if options.exportVtp
    evl.exportvtp(evl.t,evl.y,1e-3,'./vtk');
end

%% export to GMT

options.exportXyz=false;
if options.exportXyz
    evl.exportxyz(evl.t,evl.y,1e-3,'./xyz');
end

%% export to Relax

options.exportFlt=false;
if options.exportFlt
    evl.exportflt(evl.t,evl.y,1e-3,'./flt',length(evl.t));
end

%% residuals and variance reduction

% GPS data (network.dat contains a list of stations with name and coordinates)
gps.d=gps.dataVectorFromTimeSeries('./gps',time);

res=sqrt(sum((gps.d-gps.dataVectorFromModel(evl)).^2));
theta=(1-res^2/sum(gps.d.^2))*1e2;

if 0~=sum(sum(isnan(evl.y)))
    res=Inf;
end

fprintf('res=%f, theta=%2.4f%%\n',res,theta)

%% surface data simulation

% GPS data (network.dat contains a list of stations with name and coordinates)
gps=unicycle.manifold.okada85.gps('./gps/opts.dat',evl,3);
gps.d=gps.dataVectorFromTimeSeries('./gps',time);

timeIndex=ceil(length(evl.t)/2);
timeValue=evl.t(timeIndex);


options.plotGpsTimeSeries=false;
if options.plotGpsTimeSeries
    
    dm=gps.dataVectorFromModelAtTime(evl,timeValue);
    nfigure=2;
    nplots=floor(length(gps.x)/nfigure);
    nrows=3;
    ncolumns=ceil(nplots/nrows);
    vecsize=3;
    for j=0:nfigure-1
        color='b';
        figure(100+j);
        clf;
        set(gcf,'Name','GPS time series')
        
        for i=0:nplots-1
            stationId=floor(length(gps.x)/nfigure)*j+i+1;
            
            [gt,gn,ge,gd,gsn,gse,gsd]=textread(['./gps/' gps.stationName{stationId} '.dat'],...
                '%f %f %f %f %f %f %f','commentstyle','shell');
        
            subplot(nrows,ncolumns*vecsize,vecsize*i+1);
            %cla;
            hold on, box on;
            displacementComponent=1;
            plot(gt-gt(1),ge,'r');
            gps.plotTimeSeries(evl.t,stationId,displacementComponent,color);
            plot(0.5,gps.ur((stationId-1)*gps.vecsize+1,timeIndex),scale*gps.ur((stationId-1)*gps.vecsize+2,timeIndex),0,'k','linewidth',2);
            plot(0.5,scale*dm((stationId-1)*gps.vecsize+1),scale*dm((stationId-1)*gps.vecsize+2),0,'g','linewidth',2);
            set(gca,'xlim',[0 evl.t(end)]);
            title(sprintf('%s east',gps.stationName{stationId}));
            if (0==mod(i,nrows)), ylabel('displacement (m)'), end;
            if (i>=nplots-ncolumns), xlabel('time (yr)'), end;
            grid on
            set(gca, 'LooseInset', get(gca,'TightInset'),'xlim',[time(1) time(end)],'ylim',[-0.01 0.1])
            
            subplot(nrows,ncolumns*vecsize,vecsize*i+2);
            %cla;
            hold on, box on;
            displacementComponent=2;
            plot(gt-gt(1),gn,'r');
            gps.plotTimeSeries(evl.t,stationId,displacementComponent,color);
            set(gca,'xlim',[0 evl.t(end)]);
            title(sprintf('north'));
            if (i>=nplots-ncolumns), xlabel('time (yr)'), end;
            grid on
            set(gca, 'LooseInset', get(gca,'TightInset'),'xlim',[time(1) time(end)],'ylim',[-0.01 0.1])
            
            if (vecsize>2)
                subplot(nrows,ncolumns*vecsize,vecsize*i+3);
                %cla;
                hold on, box on;
                displacementComponent=3;
                plot(gt-gt(1),-gd,'r');
                gps.plotTimeSeries(evl.t,stationId,displacementComponent,color);
                title(sprintf('up'));
                if (i>=nplots-ncolumns), xlabel('time (yr)'), end;
                grid on
                set(gca, 'LooseInset', get(gca,'TightInset'),'xlim',[time(1) time(end)],'ylim',[-0.04 0.04])
            end
        end
    end
end


%% map

options.plotGpsMap=false;
if options.plotGpsMap
    scale=2e6;
    figure(2);clf;set(gcf,'Name','GPS displacement maps')
    hold on; box on
    timeIndex=ceil(length(evl.t)/2)-81;
    timeValue=evl.t(timeIndex);
    xlim=[-300 800]*1e3;
    ylim=[-400 500]*1e3;
    
    plot(gps.x(:,1),gps.x(:,2),'^');
    rcv.plotPatch()
    %rcv.plotPatch(evl.evt{1}.Fsd*(evt{1}.slip.*cosd(evt{1}.rake))); % stress change
    %rcv.plotPatch(-evl.y(1:evl.dgf:end,timeIndex)); % strike slip
    for stationId=1:length(gps.stationName)
        [gt,gn,ge,gd,gsn,gse,gsd]=textread(['./gps/' gps.stationName{stationId} '.dat'],...
        '%f %f %f %f %f %f %f','commentstyle','shell');
    
        [~,index]=min(abs(gt-gt(1)-timeValue));
        
        scatter(gps.x(stationId,1),gps.x(stationId,2),2e2,gps.ur((stationId-1)*gps.vecsize+3,timeIndex),'filled','MarkerEdgeColor','k');
        scatter(gps.x(stationId,1),gps.x(stationId,2),8e2,-gd(index),'filled','MarkerEdgeColor','k');
        
        
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*ge(index),scale*gn(index),'r','linewidth',2);
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*gps.ur((stationId-1)*gps.vecsize+1,timeIndex),scale*gps.ur((stationId-1)*gps.vecsize+2,timeIndex),'k','linewidth',2);
        
        h=plot.ellipse(scale*0.01,scale*0.01,0,gps.x(stationId,1)+scale*ge(index),gps.x(stationId,2)+scale*gn(index),'r');
        set(h,'linewidth',2)
        
        text(gps.x(stationId,1)+2e3,gps.x(stationId,2)-1e4,gps.stationName{stationId},'horizontalAlignment','right','verticalAlignment','top');
    end
    evt{1}.plotPatch()
    unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
    unicycle.plot.plot_faults(Fz_km*1e3,Fz_dim,xlim,ylim,[1 0.1 0.1]);
    unicycle.plot.plot_faults(sunda_km*1e3,sunda_dim,xlim,ylim,[0.1 0.1 0.1]);
    unicycle.plot.plot_faults(sumatra_km*1e3,sumatra_dim,xlim,ylim,[0.1 0.1 0.1]);
    xlabel('east (m)'), ylabel('north (m)')
    axis equal
    %set(gca,'xlim',xlim,'ylim',ylim,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    set(gca,'xlim',xlim,'ylim',ylim,'clim',[-1 1]*0.02);
    h=colorbar();
    ylabel(h,'vertical displacement (m)')
    %axis square
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end

%% time series plot

options.plotTimeSeries=false;
if options.plotTimeSeries
    index={1,fix(mean(find(rcv.a>0)))};
    component=2;
    colors='rgbk';
    
    figure(5);clf;set(gcf,'name','time series of patch slip')
    subplot(3,1,1);cla;
    if 0<src.N
        plot(evl.t,evl.t*src.slip(1),[':k']);
    end
    subplot(3,1,2);cla;
    subplot(3,1,3);cla;
    for k=1:length(index)
        pos = (index{k}-1)*evl.dgf;
        subplot(3,1,1);hold on
        plot(evl.t,evl.y(pos+component,:),'-o','MarkerSize',2,'Color',colors(k));
        %plot(t,y(pos+2,:),'b-+');
        xlabel('time (yr)');ylabel('fault slip (m)');box on;
        
        subplot(3,1,2);hold on
        Dt=diff(evl.t);
        v=abs(diff(evl.y(pos+component,:))./Dt)*y2s;
        plot(evl.t(2:end),log10(v),'-o','MarkerSize',2,'Color',colors(k));
        plot(evl.t,log10(mean(Dt.*v)/mean(Dt)))
        xlabel('time (yr)');ylabel('log10 of slip velocity (m/s)');box on;
        
        subplot(3,1,3);hold on
        plot(evl.t,log10(abs(evl.y(pos+3,:))),'-o','MarkerSize',2,'Color',colors(k));
        xlabel('time (yr)');ylabel('log10 of shear stress (MPa)');box on;
    end
    subplot(3,1,1);
    if ~isempty(evt)
        for p=1:3
            subplot(3,1,p);
            for k=1:length(evt)
                plot(evt{k}.t0*[1 1],get(gca,'ylim'),'k-')
            end
        end
    end
end

%% plot slip residuals

options.plotSlipResiduals=false;
if options.plotSlipResiduals
    scale=2e4;
    figure(2);
    clf;
    set(gcf,'Name','slip residuals')
    hold on; box on
    timeIndex=ceil(length(evl.t)/1);
    timeValue=evl.t(timeIndex);
    xlim=[-25 25]*1e3;
    ylim=[-25 25]*1e3;
    
    d=gps.dataVectorAtTime(timeValue);
    g=gps.dataVectorFromModelAtTime(evl,timeValue);
    r=d-g;
    m=gps.H'*pinv(gps.H*gps.H')*r;
    
    plot(gps.x(:,1),gps.x(:,2),'^');
    rcv.plotPatch(m(1:rcv.N)); % strike slip
    %rcv.plotPatch(m(1+rcv.N:end)); % dip slip
    %rcv.plotPatch(evl.y(2:evl.dgf:end,timeIndex));
    for stationId=1:length(gps.stationName)
        [gt,gn,ge,gd,gsn,gse,gsd]=textread(['./gps/' gps.stationName{stationId} '.dat'],...
        '%f %f %f %f %f %f %f','commentstyle','shell');
        
        [~,index]=min(abs(gt-gt(1)-timeValue));
        
        %scatter(gps.x(stationId,1),gps.x(stationId,2),2e2,gps.ur((stationId-1)*gps.vecsize+3,timeIndex)*1e3,'filled','MarkerEdgeColor','k');
        %scatter(gps.x(stationId,1),gps.x(stationId,2),8e2,-gd(index)*1e3,'filled','MarkerEdgeColor','k');
        
        % plot model
        %quiver(gps.x(stationId,1),gps.x(stationId,2),scale*gps.ur((stationId-1)*gps.vecsize+1,timeIndex),scale*gps.ur((stationId-1)*gps.vecsize+2,timeIndex),0,'k','linewidth',2);
        %quiver(gps.x(stationId,1),gps.x(stationId,2),scale*dm((stationId-1)*gps.vecsize+1),scale*dm((stationId-1)*gps.vecsize+2),0,'g','linewidth',2);
        
        % plot data
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*ge(index),scale*gn(index),0,'r','linewidth',2);
        h=plot.ellipse(scale*0.01,scale*0.01,0,gps.x(stationId,1)+scale*ge(index),gps.x(stationId,2)+scale*gn(index),'k');
        %quiver(gps.x(stationId,1),gps.x(stationId,2),scale*d((stationId-1)*gps.vecsize+1),scale*d((stationId-1)*gps.vecsize+2),0,'r','linewidth',2);
        %h=plot.ellipse(scale*0.01,scale*0.01,0,gps.x(stationId,1)+scale*d((stationId-1)*gps.vecsize+1),gps.x(stationId,2)+scale*d((stationId-1)*gps.vecsize+2),'r');
        set(h,'linewidth',2,'color','r');
        
        % plot residuals
        quiver(gps.x(stationId,1),gps.x(stationId,2),scale*r((stationId-1)*gps.vecsize+1),scale*r((stationId-1)*gps.vecsize+2),0,'g','linewidth',2);
        h=plot.ellipse(scale*0.01,scale*0.01,0,gps.x(stationId,1)+scale*r((stationId-1)*gps.vecsize+1),gps.x(stationId,2)+scale*r((stationId-1)*gps.vecsize+2),'g');
        set(h,'linewidth',2,'color','g');
        
        text(gps.x(stationId,1)+2e2,gps.x(stationId,2)-1e2,gps.stationName{stationId},'horizontalAlignment','left','verticalAlignment','top');
    end
    unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
    unicycle.plot.plot_faults(ca_faults_km*1e3,ca_faults_dim,xlim,ylim,[1 0.1 0.1]);
    xlabel('east (m)'), ylabel('north (m)')
    axis equal
    set(gca,'view',[-70 40],'xlim',xlim,'ylim',ylim)
    set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
    %set(gca,'clim',[-1 1]*1);
    h=colorbar();
    ylabel(h,'slip (m)')
    %axis square
    title(sprintf('surface displacements at time %f yr',evl.t(timeIndex)));
end

%% static properties 3d plot

options.plotGeometry=true;
if options.plotGeometry
    figure(6);clf;
    set(gcf,'name','source and receiver fault geometry')
    hold on
    component=1; % strike slip
    %component=2; % dip slip
    timeIndex=ceil(length(evl.t));
    evt{1}.plotPatch();
    %toplot=sqrt(evl.y(component:evl.dgf:end,timeIndex).^2+evl.y(2:evl.dgf:end,timeIndex).^2);
    toplot=-evl.y(component:evl.dgf:end,timeIndex);
    rcv.plotPatch(toplot);
    %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,timeIndex)',evl.y(2:evl.dgf:end,timeIndex)',5e6)
    %rcv.plotpatch(G*rcv.l./((rcv.b-rcv.a).*rcv.sigma)/1e3)
    %rcv.plotpatch(rcv.a-rcv.b)
    h=colorbar();
    ylabel(h,'afterslip (m)')
    
    unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
    unicycle.plot.plot_faults(ca_faults_km*1e3,ca_faults_dim,xlim,ylim,[1 0.1 0.1]);
    
    title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',t(k),t(end)));
    axis equal tight; %view([-30,20]);
    set(gca,'xlim',xlim,'ylim',ylim);
    clim=get(gca,'clim');
    set(gca,'clim',[-1 1]*max(abs(clim(:))));
    %set(gca,'clim',[-1 1]*0.3);
    box on, view([-70 40]), xlabel('x (km)'), ylabel('y (km)')
end



%% slip evolution movie

options.plotEvolution=false;
if options.plotEvolution
    % bounds
    xlim=[-300 800]*1e3;
    ylim=[-400 500]*1e3;
    sc=1e1;
    %component=1; % strike slip
    component=2; % dip slip
    
    figure(7);clf;
    for k=length(evl.t):20:length(evl.t)
        figure(7);
        subplot(2,1,1);cla;hold on
        toplot=evl.y(component:evl.dgf:end,k);
        rcv.plotPatch(toplot);
        %rcv.plotSlipVectors(evl.y(1:evl.dgf:end,k)',evl.y(2:evl.dgf:end,k)',1e1);
        
        unicycle.plot.plot_faults(coasts_km*1e3,coasts_dim,xlim,ylim,[0 0 0]);
        unicycle.plot.plot_faults(Fz_km*1e3,Fz_dim,xlim,ylim,[1 0.1 0.1]);
        unicycle.plot.plot_faults(sunda_km*1e3,sunda_dim,xlim,ylim,[0.1 0.1 0.1]);
        unicycle.plot.plot_faults(sumatra_km*1e3,sumatra_dim,xlim,ylim,[0.1 0.1 0.1]);
        
        h=colorbar();
        ylabel(h,'afterslip (m)');
        title(sprintf('afterslip (m) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-25 0]*1e3);
        %set(gca,'clim',[-1 1]*(max(abs(toplot))+1e-9));
        set(gca,'clim',[0 0.1]);
        box on, view([45 48]), xlabel('x (m)'), ylabel('y (m)')
        
        subplot(2,1,2);cla
        toplot=log10((evl.y(component:evl.dgf:end,k)-evl.y(component:evl.dgf:end,k-1))/(evl.t(k)-evl.t(k-1))*y2s);
        rcv.plotPatch(toplot);
        
        h=colorbar();
        ylabel(h,'velocity in strike direction (m/s)');
        title(sprintf('velocity in dip direction (m/s) t=%2.1e/%2.1e yr',evl.t(k),evl.t(end)));
        axis equal tight;
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-25 0]*1e3);
        set(gca,'clim',[-12 -1]);
        box on, view([45 48]), xlabel('x (m)'), ylabel('y (m)')
        
        pause(0.0125)
    end
end
