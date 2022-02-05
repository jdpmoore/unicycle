classdef Plotting
    properties
        style
        V
        downdip
        evl
        KO
        LO
        x
        dt
        Vmax
        slip
        gravity
        G33
    end
    methods
        function obj = Plotting(style,downdip,evl,xlocs,gravity,G33)
            V=repmat(evl.flt.Vo',size(evl.y,1)-1,1).*exp(evl.y(2:end,4:evl.flt.dgf:evl.flt.N*evl.flt.dgf));
            if ~exist(strcat(evl.knlpath,'KO_s.grd'),'file')
                [Ks,Kd]=unicycle.greens.computeDisplacementKernelsOkada85(evl.flt,evl.flt.earthModel.nu,xlocs,3);
                obj.KO={Ks,Kd};
                unicycle.export.grdwrite([0 1], [0 1],Ks,[evl.knlpath,'KO_s.grd'])
                unicycle.export.grdwrite([0 1], [0 1],Kd,[evl.knlpath,'KO_d.grd'])
            else
                [~,~,Ks]=unicycle.export.grdread([evl.knlpath,'KO_s.grd']);
                [~,~,Kd]=unicycle.export.grdread([evl.knlpath,'KO_d.grd']);
                obj.KO={Ks,Kd};
            end
            if evl.shz.N > 0
                if ~exist(strcat(evl.knlpath,'LO_11.grd'),'file')
                    [L11,L12,L13,L22,L23,L33]=computeDisplacementKernelsVerticalShearZone(evl.shz,evl.shz.earthModel.nu,xlocs,3);
                    obj.LO={L11,L12,L13,L22,L23,L33};
                    unicycle.export.grdwrite([0 1], [0 1],L11,[evl.knlpath,'LO_11.grd'])
                    unicycle.export.grdwrite([0 1], [0 1],L12,[evl.knlpath,'LO_12.grd'])
                    unicycle.export.grdwrite([0 1], [0 1],L13,[evl.knlpath,'LO_13.grd'])
                    unicycle.export.grdwrite([0 1], [0 1],L22,[evl.knlpath,'LO_22.grd'])
                    unicycle.export.grdwrite([0 1], [0 1],L23,[evl.knlpath,'LO_23.grd'])
                    unicycle.export.grdwrite([0 1], [0 1],L33,[evl.knlpath,'LO_33.grd'])
                else
                    [~,~,L11]=unicycle.export.grdread([evl.knlpath,'LO_11.grd']);
                    [~,~,L12]=unicycle.export.grdread([evl.knlpath,'LO_12.grd']);
                    [~,~,L13]=unicycle.export.grdread([evl.knlpath,'LO_13.grd']);
                    [~,~,L22]=unicycle.export.grdread([evl.knlpath,'LO_22.grd']);
                    [~,~,L23]=unicycle.export.grdread([evl.knlpath,'LO_23.grd']);
                    [~,~,L33]=unicycle.export.grdread([evl.knlpath,'LO_33.grd']);
                    obj.LO={L11,L12,L13,L22,L23,L33};
                end
            else
                obj.LO={0,0,0,0,0,0};
            end
            obj.style = style;
            obj.V = V;
            obj.Vmax=max(V,[],2);
            obj.downdip = downdip;
            obj.evl = evl;
            obj.x = xlocs;
            obj.gravity = gravity;
            obj.G33 = G33;
            %Yp=zeros(length(evl.t)-1,size(evl.y,2));
            dt=zeros(length(evl.t)-1,1);
            for k=1:length(evl.t)-1
                %Yp(k,:)=(evl.y(k+1,:)-evl.y(k,:))/(evl.t(k+1)-evl.t(k));
                dt(k)=(evl.t(k+1)-evl.t(k));
            end
            obj.dt=dt;
            obj.slip = evl.y(2:end,1:evl.flt.dgf:evl.flt.N*evl.flt.dgf)';
        end
        %
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % %                    F I G U R E S                     %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %  %
        %
        % if style == 1
        %     figure(1011);clf;set(gcf,'Color','White','name','Time Evolution')
        %     f11a = subplot(2,1,1);cla;
        %     toplot=Ks*evl.y(:,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %surface displacements
        %     toplot=toplot(2:3:end,:);
        %     pcolor(evl.t/3.15e7,xloc/1e3,toplot), shading flat
        %     set(gca,'YDir','reverse');
        %     h=colorbar('Location','NorthOutside');
        %     %caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
        %     caxis(minmax(toplot));
        %     colormap(f11a,jet);
        %     title(h,'Surface deformation (m)')
        %     xlabel('Time (yr)')
        %     ylabel('Distance from fault (E-W km)');
        %     f11b = subplot(2,1,2);cla;
        %     pcolor(1:length(evl.t),xloc/1e3,toplot), shading flat
        %     set(gca,'YDir','reverse');
        %     h=colorbar('Location','NorthOutside');
        %     caxis(minmax(toplot));
        %     colormap(f11b,jet);
        %     title(h,'Surface deformation (m)')
        %     xlabel('Time steps')
        %     ylabel('Distance from fault (E-W km)');
        % end
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % %               Time Tracking of Points                %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        % if style == 1
        % figure(1012);clf;set(gcf,'Color','White','name','Time evolution of tracked points');gca; hold on;box on;grid on;
        % end
        % if style == 2
        %     toplot = Kd*evl.y(:,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %     toplotE = toplot(1:3:end,:);
        %     toplotU = toplot(3:3:end,:);
        %     figure(1012);clf;set(gcf,'Color','White','name','Time evolution of tracked points');subplot(2,1,1); gca; hold on;box on;grid on;
        %     plot(evl.t./year,toplotE(600,:))
        %     subplot(2,1,2); gca; hold on;box on;grid on;
        %     plot(evl.t./year,toplotU(600,:))
        % end
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % %          CONTOURS OF SURFACE DISPLACEMENTS           %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        % if true
        %     if style == 1
        %         figure(1009);clf;set(gcf,'Color','White','name','Surface deformation');gca; hold on;box on;grid on;
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Ks*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(2:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/52);
        %             toplots2=Ks*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(2:3:end,:);
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),25*year); % aseismic slip
        %             toplots3=Ks*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(2:3:end,:);
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Along-strike displacements (m)');
        %         title('Along-strike surface displacements');
        % %        lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (1s)','Postseismic (1m)','Interseismic (25yr)'},'FontSize',18);
        %  %       legend('boxoff')
        % %        title(lgd,'Contour spacing')
        %     end
        %     if style == 2
        %         figure(1009);clf;set(gcf,'Color','White','name','Surface deformation');subplot(2,1,1);gca; hold on;box on;grid on;
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Kd*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(1:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/12);
        %             toplots2=Kd*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(1:3:end,:);
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),5*year); % aseismic slip
        %             toplots3=Kd*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(1:3:end,:);
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Horizontal displacement (m)');
        %         title('Convergence');
        % %        lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (100ms)','Postseismic (1m)','Interseismic (5yr)'},'FontSize',18);
        % %        legend('boxoff')
        % %        title(lgd,'Contour spacing')
        %
        %         subplot(2,1,2);gca; hold on;box on;grid on;
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Kd*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(3:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/12);
        %             toplots2=Kd*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(3:3:end,:);
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),5*year); % aseismic slip
        %             toplots3=Kd*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(3:3:end,:);
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Vertical displacement (m)');
        %         title('Uplift');
        % %         lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (100ms)','Postseismic (1m)','Interseismic (5yr)'},'FontSize',18);
        % %         legend('boxoff')
        % %         title(lgd,'Contour spacing')
        %     end
        % end
        %
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % % CONTOURS OF SURFACE DISPLACEMENTS GRAVITY CORRECTED  %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        % if true
        %     if style == 1
        %         figure(1011);clf;set(gcf,'Color','White','name','Surface deformation');gca; hold on;box on;grid on;
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Ks*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(2:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/52);
        %             toplots2=Ks*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(2:3:end,:);
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),25*year); % aseismic slip
        %             toplots3=Ks*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(2:3:end,:);
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Along-strike displacements (m)');
        %         title('Along-strike surface displacements');
        % %         lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (1s)','Postseismic (1m)','Interseismic (25yr)'},'FontSize',18);
        % %         legend('boxoff')
        % %         title(lgd,'Contour spacing')
        %     end
        %     if style == 2
        %         figure(1011);clf;set(gcf,'Color','White','name','Surface deformation');subplot(2,1,1);gca; hold on;box on;grid on;
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Kd*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(1:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/12);
        %             toplots2=Kd*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(1:3:end,:);
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),5*year); % aseismic slip
        %             toplots3=Kd*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(1:3:end,:);
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Horizontal displacement (m)');
        %         title('Convergence');
        % %         lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (100ms)','Postseismic (1m)','Interseismic (5yr)'},'FontSize',18);
        % %         legend('boxoff')
        % %         title(lgd,'Contour spacing')
        %
        %         subplot(2,1,2);gca; hold on;box on;grid on;
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.1);
        %             toplots1=Kd*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %seismic displacements
        %             toplots1=toplots1(3:3:end,:);
        %             toplots1=toplots1+g*G33*toplots1;
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/12);
        %             toplots2=Kd*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots2=toplots2(3:3:end,:);
        %             toplots2=toplots2+g*G33*toplots2;
        %             p2=plot(xloc/1e3,toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),5*year); % aseismic slip
        %             toplots3=Kd*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=toplots3(3:3:end,:);
        %             toplots3=toplots3+g*G33*toplots3./1e5;
        %             p3=plot(xloc/1e3,toplots3,'b','LineWidth',1.5);
        %         end
        %         axis tight
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Vertical displacement (m)');
        %         title('Uplift');
        % %         lgd = legend([p1(1) p2(1) p3(1)],{'Coseismic (100ms)','Postseismic (1m)','Interseismic (5yr)'},'FontSize',18);
        % %         legend('boxoff')
        % %         title(lgd,'Contour spacing')
        %     end
        % end
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % %            CONTOURS OF SURFACE VELOCITIES            %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        % if true
        %     if style == 1
        %         figure(1010);clf;set(gcf,'Color','White','name','Surface deformation');subplot(3,1,1);gca; hold on;box on;grid on;
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.01);
        %             deltat=evl.t(V1index(pos))-evl.t(V1index(pos)-1);
        %             toplots1u=Ks*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots1l=Ks*evl.y(V1index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots1=(toplots1u-toplots1l)./deltat';
        %             toplots1=toplots1(2:3:end,:);
        %             p1=plot(xloc/1e3,toplots1,'r','LineWidth',1.0);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Along-strike velocities (m/s)');
        %             title('Along-strike co-seismic surface velocities');
        %         end
        %         if(~isempty(V2index))
        %             subplot(3,1,2);gca; hold on;box on;grid on;
        %             pos = getPeriodicIndex(evl.t(V2index),year/365);
        %             deltat=evl.t(V2index(pos))-evl.t(V2index(pos)-1);
        %             toplots2u=Ks*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots2l=Ks*evl.y(V2index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots2=(toplots2u-toplots2l)./deltat';
        %             toplots2=toplots2(2:3:end,:);
        %             p2=plot(xloc/1e3,year*toplots2,'color',[0 0.7 0],'LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Along-strike velocities (m/yr)');
        %             title('Along-strike post-seismic surface velocities');
        %         end
        %         if(~isempty(V3index))
        %             subplot(3,1,3);gca; hold on;box on;grid on;
        %             pos = getPeriodicIndex(evl.t(V3index),1*year); % aseismic slip
        %             if pos(1) == 1
        %                 pos = pos(2:end);
        %             end
        %             deltat=evl.t(V3index(pos))-evl.t(V3index(pos)-1);
        %             toplots3u=Ks*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3l=Ks*evl.y(V3index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=(toplots3u-toplots3l)./deltat';
        %             toplots3=toplots3(2:3:end,:);
        %             p3=plot(xloc/1e3,1e3*year*toplots3,'b','LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Along-strike velocities (mm/yr)');
        %             title('Along-strike inter-seismic surface velocities');
        %         end
        %     end
        %     if style == 2
        %         figure(1010);clf;set(gcf,'Color','White','name','Surface deformation');
        %         Vmaxl = log10(Vmax);
        %         V1index= find(Vmaxl>=-2); % seismic velocity
        %         V2index= find(ceil(log10(max(flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
        %         V3index= find(Vmaxl<ceil(log10(max(flt.Vpl)))); % aseismic velocity
        %         if (~isempty(V1index))
        %             pos = getPeriodicIndex(evl.t(V1index),0.01);
        %             deltat=evl.t(V1index(pos))-evl.t(V1index(pos)-1);
        %             toplots1u=Kd*evl.y(V1index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots1l=Kd*evl.y(V1index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots1=(toplots1u-toplots1l)./deltat';
        %             subplot(3,2,1);gca; hold on;box on;grid on;
        %             toplots1H=toplots1(1:3:end,:);
        %             p1=plot(xloc/1e3,toplots1H,'r','LineWidth',1.0);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Horizontal velocities (m/s)');
        %             title('Co-seismic surface horizontal velocities');
        %             subplot(3,2,2);gca; hold on;box on;grid on;
        %             toplots1V=toplots1(3:3:end,:);
        %             p1=plot(xloc/1e3,toplots1V,'r','LineWidth',1.0);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Vertical velocities (m/s)');
        %             title('Co-seismic surface uplift velocities');
        %         end
        %         if(~isempty(V2index))
        %             pos = getPeriodicIndex(evl.t(V2index),year/365);
        %             deltat=evl.t(V2index(pos))-evl.t(V2index(pos)-1);
        %             toplots2u=Kd*evl.y(V2index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots2l=Kd*evl.y(V2index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplots2=(toplots2u-toplots2l)./deltat';
        %             subplot(3,2,3);gca; hold on;box on;grid on;
        %             toplots2H=toplots2(1:3:end,:);
        %             p2=plot(xloc/1e3,year*toplots2H,'color',[0 0.7 0],'LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Horizontal velocities (m/yr)');
        %             title('Post-seismic surface horizontal velocities');
        %             subplot(3,2,4);gca; hold on;box on;grid on;
        %             toplots2V=toplots2(3:3:end,:);
        %             p2=plot(xloc/1e3,year*toplots2V,'color',[0 0.7 0],'LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Vertical velocities (m/yr)');
        %             title('Post-seismic surface uplift velocities');
        %         end
        %         if(~isempty(V3index))
        %             pos = getPeriodicIndex(evl.t(V3index),1*year); % aseismic slip
        %             if pos(1) == 1
        %                 pos = pos(2:end);
        %             end
        %             deltat=evl.t(V3index(pos))-evl.t(V3index(pos)-1);
        %             toplots3u=Kd*evl.y(V3index(pos),1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3l=Kd*evl.y(V3index(pos)-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))'; %postseismic displacements
        %             toplots3=(toplots3u-toplots3l)./deltat';
        %             subplot(3,2,5);gca; hold on;box on;grid on;
        %             toplots3H=toplots3(1:3:end,:);
        %             p3=plot(xloc/1e3,1e3*year*toplots3H,'b','LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Horizontal velocities (mm/yr)');
        %             title('Inter-seismic surface horizontal velocities');
        %             subplot(3,2,6);gca; hold on;box on;grid on;
        %             toplots3V=toplots3(3:3:end,:);
        %             p3=plot(xloc/1e3,1e3*year*toplots3V,'b','LineWidth',1.5);
        %             axis tight
        %             xlabel('Distance from fault (E-W km)');
        %             ylabel('Vertical velocities (mm/yr)');
        %             title('Inter-seismic surface uplift velocities');
        %         end
        %     end
        % end
        % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % %                                                      %
        % %       SURFACE DISPLACEMENTS AND VELOCITIES V1        %
        % %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        % colours='mbrkcy';
        % if true
        %     numcurves=5;
        %     if style==1
        %         figure(1007);clf;set(gcf,'Color','White','name','Surface deformation');gca; hold on;box on;grid on;
        %         toplotV=zeros(numcurves,length(xloc));
        %         for k=1:numcurves
        %             tindex=round(k*(length(evl.t)/numcurves));
        %             toplota=Ks*evl.y(tindex,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplotVE=L12*evl.y(tindex,(evl.flt.N*evl.flt.dgf+1):evl.shz.dgf:end)'+L13*evl.y(tindex,(evl.flt.N*evl.flt.dgf+2):evl.shz.dgf:end)';
        %             toplotaV=(toplota-Ks*evl.y(tindex-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))')./(evl.t(tindex)-evl.t(tindex-1));
        %             toplota=toplota(2:3:end);
        %             toplotVE=2e1*toplotVE(2:3:end);
        %             toplotV(k,:)=toplotaV(2:3:end);
        %             plot(xloc./1e3,toplota,colours(k),'LineWidth',2.0)
        %             plot(xloc./1e3,toplotVE,[colours(k) '--'],'LineWidth',2.0)
        %             plot(xloc./1e3,toplota+toplotVE,[colours(k) '-.'],'LineWidth',2.0)
        %             legendnames{k}=num2str(round(evl.t(tindex)./year));
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Along-strike deformation (m)');
        %         axis tight;
        %         title('Surface deformation');
        %         lgd = legend(legendnames,'FontSize',14);
        %         legend('boxoff')
        %         title(lgd,'Time (yrs)')
        %         figure(1008);clf;set(gcf,'Color','White','name','Surface velocities');gca; hold on;box on;grid on;
        %         for k=1:numcurves
        %             plot(xloc./1e3,toplotV(k,:),'LineWidth',2.0)
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Along-strike velocities (m/s)');
        %         title('Surface velocities');
        %         lgd = legend(legendnames,'FontSize',14);
        %         legend('boxoff');axis tight;
        %         title(lgd,'Time (yrs)')
        %     end
        %     if style==2
        %         figure(1007);clf;set(gcf,'Color','White','name','Surface deformation');subplot(2,1,1);gca; hold on;box on;grid on;
        %         toplotVc=zeros(numcurves,length(xloc));
        %         toplotVu=zeros(numcurves,length(xloc));
        %         for k=1:numcurves
        %             tindex=round(k*(length(evl.t)/numcurves));
        %             toplota=Kd*evl.y(tindex,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplotaV=(toplota-Kd*evl.y(tindex-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))')./(evl.t(tindex)-evl.t(tindex-1));
        %             toplotac=toplota(1:3:end);
        %             toplotVc(k,:)=toplotaV(1:3:end);
        %             plot(xloc./1e3,toplotac,colors(k),'LineWidth',2.0)
        %             legendnames{k}=num2str(round(evl.t(tindex)./year));
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Horizontal deformation (m)');
        %         title('Convergence');
        %         lgd = legend(legendnames,'FontSize',10);
        %         legend('boxoff')
        %         title(lgd,'Time (yrs)')
        %         subplot(2,1,2);gca; hold on;box on;grid on;
        %         for k=1:numcurves
        %             tindex=round(k*(length(evl.t)/numcurves));
        %             toplota=Kd*evl.y(tindex,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))';
        %             toplotaV=(toplota-Kd*evl.y(tindex-1,1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf))')./(evl.t(tindex)-evl.t(tindex-1));
        %             toplotau=toplota(3:3:end);
        %             toplotGrav=toplotau+g*G33*toplotau;
        %             %toplotGrav=G33*toplotau./1e5;
        %             toplotVu(k,:)=toplotaV(3:3:end);
        %             plot(xloc./1e3,toplotau,colors(k),'LineWidth',2.0)
        %             plot(xloc./1e3,toplotGrav,['--' colors(k)],'LineWidth',2.0)
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Vertical deformation (m)');
        %         title('Uplift');
        %         lgd = legend(legendnames,'FontSize',10);
        %         legend('boxoff')
        %         title(lgd,'Time (yrs)')
        %         figure(1008);clf;set(gcf,'Color','White','name','Surface velocities');subplot(2,1,1);gca; hold on;box on;grid on;
        %         for k=1:numcurves
        %             plot(xloc./1e3,toplotVc(k,:),'LineWidth',2.0)
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Convergence velocity (m/s)');
        %         title('Convergence velocity');
        %         lgd = legend(legendnames,'FontSize',10);
        %         legend('boxoff')
        %         title(lgd,'Time (yrs)')
        %         subplot(2,1,2);gca; hold on;box on;grid on;
        %         for k=1:numcurves
        %             plot(xloc./1e3,toplotVu(k,:),'LineWidth',2.0)
        %         end
        %         xlabel('Distance from fault (E-W km)');
        %         ylabel('Uplift velocity (m/s)');
        %         title('Uplift velocities');
        %         lgd = legend(legendnames,'FontSize',10);
        %         legend('boxoff')
        %         title(lgd,'Time (yrs)')
        %     end
        % end
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                      %
        %       SURFACE DISPLACEMENTS AND VELOCITIES V1        %
        %                                                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        
        function surfaceDisplacement(obj,fignum,numcurves)
            year = 60*60*24*365;g=9.81;
            colours='mbrkcy';
            if obj.style==1
                figure(fignum(1));clf;set(gcf,'Color','White','name','Surface deformation');gca; hold on;box on;grid on;
                toplotV=zeros(numcurves,length(obj.x));
                for k=1:numcurves
                    tindex=round(k*(length(obj.evl.t)/numcurves));
                    toplota=obj.KO{1}*obj.evl.y(tindex,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))';
                    toplotVE=obj.LO{2}*obj.evl.y(tindex,(obj.evl.flt.N*obj.evl.flt.dgf+1):obj.evl.shz.dgf:end)'+obj.LO{3}*obj.evl.y(tindex,(obj.evl.flt.N*obj.evl.flt.dgf+2):obj.evl.shz.dgf:end)';
                    toplotaV=(toplota-obj.KO{1}*obj.evl.y(tindex-1,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))')./(obj.evl.t(tindex)-obj.evl.t(tindex-1));
                    toplota=toplota(2:3:end);
                    toplotVE=2e1*toplotVE(2:3:end);
                    toplotV(k,:)=toplotaV(2:3:end);
                    plot(obj.x(:,1)./1e3,toplota,colours(k),'LineWidth',2.0)
                    if size(toplotVE,1) > 0
                        plot(obj.x(:,1)./1e3,toplotVE,[colours(k) '--'],'LineWidth',2.0)
                        plot(obj.x(:,1)./1e3,toplota+toplotVE,[colours(k) '-.'],'LineWidth',2.0)
                    end
                    legendnames{k}=num2str(round(obj.evl.t(tindex)./year));
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Along-strike deformation (m)');
                axis tight;
                title('Surface deformation');
                lgd = legend(legendnames,'FontSize',14);
                legend('boxoff')
                title(lgd,'Time (yrs)')
                figure(fignum(2));clf;set(gcf,'Color','White','name','Surface velocities');gca; hold on;box on;grid on;
                for k=1:numcurves
                    plot(obj.x(:,1)./1e3,toplotV(k,:),'LineWidth',2.0)
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Along-strike velocities (m/s)');
                title('Surface velocities');
                lgd = legend(legendnames,'FontSize',14);
                legend('boxoff');axis tight;
                title(lgd,'Time (yrs)')
            end
            if obj.style==2
                figure(fignum(1));clf;set(gcf,'Color','White','name','Surface deformation');subplot(2,1,1);gca; hold on;box on;grid on;
                toplotVc=zeros(numcurves,length(obj.x));
                toplotVu=zeros(numcurves,length(obj.x));
                for k=1:numcurves
                    tindex=round(k*(length(obj.evl.t)/numcurves));
                    toplota=obj.KO{2}*obj.evl.y(tindex,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))';
                    toplotaV=(toplota-obj.KO{2}*obj.evl.y(tindex-1,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))')./(obj.evl.t(tindex)-obj.evl.t(tindex-1));
                    toplotac=toplota(1:3:end);
                    toplotau=toplota(3:3:end);
                    toplotVc(k,:)=toplotaV(1:3:end);
                    plot(obj.x(:,1)./1e3,toplotac,colours(k),'LineWidth',2.0)
                    if obj.gravity
                        toplotGrav=toplotac+g*obj.G33(1:3:end,:)*toplotau;
                        %toplotGrav=g*obj.G33(3:3:end,:)*toplotau;
                        plot(obj.x(:,1)./1e3,toplotGrav,['--' colours(k)],'LineWidth',2.0)
                    end
                    legendnames{k}=num2str(round(obj.evl.t(tindex)./year));
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Horizontal deformation (m)');
                title('Convergence');
                lgd = legend(legendnames,'FontSize',10);
                legend('boxoff')
                title(lgd,'Time (yrs)')
                subplot(2,1,2);gca; hold on;box on;grid on;
                for k=1:numcurves
                    tindex=round(k*(length(obj.evl.t)/numcurves));
                    toplota=obj.KO{2}*obj.evl.y(tindex,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))';
                    toplotaV=(toplota-obj.KO{2}*obj.evl.y(tindex-1,1:obj.evl.flt.dgf:(obj.evl.flt.N*obj.evl.flt.dgf))')./(obj.evl.t(tindex)-obj.evl.t(tindex-1));
                    toplotau=toplota(3:3:end);
                    toplotVu(k,:)=toplotaV(3:3:end);
                    plot(obj.x(:,1)./1e3,toplotau,colours(k),'LineWidth',2.0)
                    if obj.gravity
                        toplotGrav=toplotau+g*obj.G33(3:3:end,:)*toplotau;
                        %toplotGrav=g*obj.G33(3:3:end,:)*toplotau;
                        plot(obj.x(:,1)./1e3,toplotGrav,['--' colours(k)],'LineWidth',2.0)
                    end
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Vertical deformation (m)');
                title('Uplift');
                lgd = legend(legendnames,'FontSize',10);
                legend('boxoff')
                title(lgd,'Time (yrs)')
                figure(fignum(2));clf;set(gcf,'Color','White','name','Surface velocities');subplot(2,1,1);gca; hold on;box on;grid on;
                for k=1:numcurves
                    plot(obj.x(:,1)./1e3,toplotVc(k,:),'LineWidth',2.0)
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Convergence velocity (m/s)');
                title('Convergence velocity');
                lgd = legend(legendnames,'FontSize',10);
                legend('boxoff')
                title(lgd,'Time (yrs)')
                subplot(2,1,2);gca; hold on;box on;grid on;
                for k=1:numcurves
                    plot(obj.x(:,1)./1e3,toplotVu(k,:),'LineWidth',2.0)
                end
                xlabel('Distance from fault (E-W km)');
                ylabel('Uplift velocity (m/s)');
                title('Uplift velocities');
                lgd = legend(legendnames,'FontSize',10);
                legend('boxoff')
                title(lgd,'Time (yrs)')
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                Slip Contours                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function slipContours(obj,fignum,periods,toSave)
            year = 60*60*24*365;
            f5 = figure(fignum);clf;set(gcf,'Color','White','name','Accumulated Slip')
            Vmaxl = log10(obj.Vmax);
            cla;hold on;
            V1index= find(Vmaxl>=-2); % seismic velocity
            V2index= find(ceil(log10(max(obj.evl.flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
            V3index= find(Vmaxl<ceil(log10(max(obj.evl.flt.Vpl)))); % aseismic velocity
            leg = [];
            leglab = {};
            if (~isempty(V1index))
                pos = getPeriodicIndex(obj.evl.t(V1index),periods(1));
                s1= obj.slip(:,V1index(pos)); % seismic slip
                p1=plot(s1,obj.downdip,'r','LineWidth',1.0);
                leg = [leg p1(1)];
                leglab = [leglab ['Coseismic (' num2str(periods(1)) ' seconds)']];
            end
            set(gca,'ydir','reverse')
            if(~isempty(V2index))
                pos = getPeriodicIndex(obj.evl.t(V2index),periods(2));
                s2= obj.slip(:,V2index(pos)); % afterslip
                p2=plot(s2,obj.downdip,'color',[0 0.7 0],'LineWidth',1.5);
                leg = [leg p2(1)];
                leglab = [leglab ['Postseismic (' num2str(12*periods(2)/year) ' months)']];
            end
            if(~isempty(V3index))
                pos = getPeriodicIndex(obj.evl.t(V3index),periods(3)); % aseismic slip
                s3= obj.slip(:,V3index(pos));
                p3=plot(s3,obj.downdip,'b','LineWidth',1.0);
                leg = [leg p3(1)];
                leglab = [leglab ['Interseismic (' num2str(periods(3)/year) ' years)']];
            end
            box on;grid on;axis tight;
            xlabel('Accumulated slip (m)')
            ylabel('Down-dip distance (km)');
            lgd = legend(leg,leglab,'FontSize',18,'Location','northoutside');
            legend('boxoff')
            title(lgd,'Contour spacing')
            % export slip contours for plotting in Mathematica
            if toSave==true
                downdip=obj.downdip;
                save(['./faultslipcontours.mat'],'downdip','s1','s2','s3')
            end
        end
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                Overprinting Slip Contours                      %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function slipContoursOver(obj,periods)
            year = 60*60*24*365;
            Vmaxl = log10(obj.Vmax);
            hold on;
            V1index= find(Vmaxl>=-2); % seismic velocity
            V2index= find(ceil(log10(max(obj.evl.flt.Vpl)))<=Vmaxl & Vmaxl<-2); %
            V3index= find(Vmaxl<ceil(log10(max(obj.evl.flt.Vpl)))); % aseismic velocity
            if (~isempty(V1index))
                pos = getPeriodicIndex(obj.evl.t(V1index),periods(1));
                s1= obj.slip(:,V1index(pos)); % seismic slip
                p1=plot(s1,obj.downdip,'r','LineWidth',1.0);
            end
            if(~isempty(V2index))
                pos = getPeriodicIndex(obj.evl.t(V2index),periods(2));
                s2= obj.slip(:,V2index(pos)); % afterslip
                p2=plot(s2,obj.downdip,'color',[0 0.7 0],'LineWidth',1.5);
            end
            if(~isempty(V3index))
                pos = getPeriodicIndex(obj.evl.t(V3index),periods(3)); % aseismic slip
                s3= obj.slip(:,V3index(pos));
                p3=plot(s3,obj.downdip,'b','LineWidth',1.0);
            end
        end
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                   Fault Slip Rate                   %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function velocityTime(obj,fignum,stressing)
            year = 60*60*24*365;
            f3 = figure(fignum);clf;set(gcf,'Color','White','name','Time Evolution')
            f3a = subplot(2,1,1);cla;
            pcolor(obj.evl.t(1:end-1)/year,obj.downdip,log10(obj.V')), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            %caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
            caxis([-12 1]);
            colormap(f3a,parula);
            title(h,'Log10 of Slip Rate (m/s)')
            xlabel('Time (yr)')
            ylabel('Down-dip distance (km)');
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
            %                Variable Stressing                   %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            f3b = subplot(2,1,2);cla;
            stress=obj.evl.flt.sigma+stressing(obj.evl.t');
            pcolor(obj.evl.t/year,obj.downdip,stress), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f3b,flipud(hot));
            title(h,'Fault normal stress (MPa)')
            xlabel('Time (yr)')
            ylabel('Down-dip distance (km)');
        end
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                   Erics Awesome new figure                   %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function velocityEric(obj,fignum,nxpts)
            year = 60*60*24*365;
            f3 = figure(fignum);clf;set(gcf,'Color','White','name','Time Evolution')
            cla;
            ydata=repmat(obj.downdip,[1,size(obj.evl.t(1:end-1))]);
            xdata=obj.slip;
            zdata=log10(obj.V');
            xsample=linspace(min(min(obj.slip)),max(max(obj.slip)),nxpts);
            test1=griddata(xdata(:),ydata(:),zdata(:),xsample,obj.downdip);
            pcolor(xsample,obj.downdip,test1), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            %caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
            caxis([-12 1]);
            colormap(f3,parula);
            title(h,'Log10 of Slip Rate (m/s)')
            xlabel('Accumulated slip (m)')
            ylabel('Down-dip distance (km)');
        end
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %               Size of Time Steps                    %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function velocityTimestep(obj,fignum,stressing)
            year = 60*60*24*365;
            f4 = figure(fignum);clf;set(gcf,'Color','White','name','Time-Step Evolution')
            f4a=subplot(3,1,1);cla;hold on
            plot(1:length(obj.evl.t)-1,log10(obj.dt),'b','LineWidth',1.0)
            xlabel('Time Steps')
            ylabel('log10(dt)')
            day = 24*60*60;
            yticks([-3 0 log10(60*60) log10(day) log10(7*day) log10(30*day) log10(365*day)]);
            yticklabels({'milisecond','second','hour','day','week','month','year'})
            f4a.YGrid = 'on';
            axis tight
            set(gca,'ylim',[min(log10(obj.dt)) log10(year)]);
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %      Slip Rate as a function of Time Step           %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            f4b = subplot(3,1,2);cla;
            pcolor(1:length(obj.evl.t)-1,obj.downdip,log10(obj.V')), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            %caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
            caxis([-12 1]);
            colormap(f4b,parula);
            title(h,'Log10 of Slip Rate (m/s)')
            xlabel('Time Steps')
            ylabel('Down-dip distance (km)');
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %      Variable Stressing as a function of time step  %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            f4c = subplot(3,1,3);cla;
            stress=obj.evl.flt.sigma+stressing(obj.evl.t');
            pcolor(1:length(obj.evl.t),obj.downdip,stress), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f4c,flipud(hot));
            title(h,'Fault normal stress (MPa)')
            xlabel('Time Steps')
            ylabel('Down-dip distance (km)');
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %  Export velocities to grd file for plotting in GMT %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportGMT(obj,fname)
            minmax=@(x) [min(x(:)),max(x(:))]; % Minimum and maximum function
            year = 60*60*24*365;
            tmax=obj.evl.t(end);
            tyrs=obj.evl.t(2:end)./year;
            tmaxyr=tmax/year;
            [tt,xx3]=meshgrid(tyrs,obj.downdip);
            toExport=interp2(tt,xx3,log10(obj.V'),1:0.01:tmaxyr,obj.downdip);
            unicycle.export.grdwrite([0 tmaxyr],minmax(obj.downdip),toExport,[fname '_slip.grd']);
            toExport=log10(obj.V');
            unicycle.export.grdwrite([0 length(obj.evl.t)],minmax(obj.downdip),toExport,[fname '_steps.grd']);
            unicycle.export.grdwrite([0 length(obj.evl.t)],minmax(obj.downdip),[log10(obj.dt) log10(obj.dt)]',[fname '_dt.grd']);
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                 Slip deficit plot                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function slipDeficit(obj,fignum)
            figure(fignum);clf;set(gcf,'Color','White','name','Slip deficit');gca; hold on;box on;grid on;
            endsliprate=(obj.V(end,:))/max(obj.V(end,:));
            plot(obj.downdip,endsliprate,'b','LineWidth',2.0)
            xlabel('Down-dip distance (km)');
            ylabel('Slip-rate (V/V_{pl})');
            title('Slip deficit');
            set(gca,'ylim',[0 1]);
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %          Viscoelastic deformation                   %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function viscoelastic(obj,fignum,time)
            year = 60*60*24*365;
            tindex=getClosestIndex(obj.evl.t./year, time);
            figure(fignum);clf;set(gcf,'Color','White','name','Viscous Strain')
            e12=obj.evl.y(tindex,obj.evl.flt.N*obj.evl.flt.dgf+1:obj.evl.shz.dgf:end);
            e13=obj.evl.y(tindex,obj.evl.flt.N*obj.evl.flt.dgf+2:obj.evl.shz.dgf:end);
            f3a = subplot(2,1,1);cla;
            x1u=unique(obj.evl.shz.xc(:,1))./1e3;
            x3u=-unique(obj.evl.shz.xc(:,3))./1e3;
            x3u=flip(x3u);
            len1=length(x1u);
            len2=length(x3u);
            pcolor(x1u,x3u,1e6.*reshape(e12,[len1 len2])'), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f3a,parula);
            title(h,['Micro-Strain (e_{12}) at ' num2str(round(obj.evl.t(tindex)./year)) 'yrs'])
            xlabel('Horizontal Location (km)')
            ylabel('Depth (km)');
            f3b = subplot(2,1,2);cla;
            pcolor(x1u,x3u,1e6.*reshape(e13,[len1 len2])'), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f3b,parula);
            title(h,['Micro-Strain (e_{13}) at ' num2str(round(obj.evl.t(tindex)./year)) 'yrs'])
            xlabel('Horizontal Location (km)')
            ylabel('Depth (km)');
        end
        
        function viscoelasticTime(obj,fignum,X)
            x1u=unique(obj.evl.shz.xc(:,1))./1e3;
            x3u=-unique(obj.evl.shz.xc(:,3))./1e3;
            x3u=flip(x3u);
            len1=length(x1u);
            len2=length(x3u);
            Xindex=getClosestIndex(x1u, X);
            year = 60*60*24*365;
            figure(fignum);clf;set(gcf,'Color','White','name','Viscous Strain')
            e12t=obj.evl.y(:,obj.evl.flt.N*obj.evl.flt.dgf+1:obj.evl.shz.dgf:end);
            e12tb=e12t(:,Xindex:len1:end);
            e13t=obj.evl.y(:,obj.evl.flt.N*obj.evl.flt.dgf+2:obj.evl.shz.dgf:end);
            e13tb=e13t(:,Xindex:len1:end);
            f3a = subplot(2,1,1);cla;
            pcolor(obj.evl.t/year,x3u,1e6.*e12tb'), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f3a,parula);
            title(h,['Micro-Strain (e_{12}) at ' num2str(round(obj.evl.shz.xc(Xindex,1)./1e3)) 'km'])
            xlabel('Time (yrs)')
            ylabel('Depth (km)');
            f3a = subplot(2,1,2);cla;
            pcolor(obj.evl.t/year,x3u,1e6.*e13tb'), shading flat
            set(gca,'YDir','reverse');
            h=colorbar('Location','NorthOutside');
            colormap(f3a,parula);
            title(h,['Micro-Strain (e_{13}) at ' num2str(round(obj.evl.shz.xc(Xindex,1)./1e3)) 'km'])
            xlabel('Time (yrs)')
            ylabel('Depth (km)');
        end
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                    CITATATIONS                     %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        %plotter.citations([viscoelastic gravity pinning tstresses])
        function citations(~,opts)
            fprintf('Required citations:\n------\n');
            fprintf('Sato D and Moore J D P.\nDisplacements and stress associated with localised and distributed inelastic deformation with piecewise-constant elastic variations\nGJI (2022)\n');
            fprintf('------\nOkada Y.\nInternal deformation due to shear and tensile faults in a half-space\nBSSA (1992)\n');
            fprintf('------\nJames D P Moore, Sylvain Barbot, Lujia Feng, Yu Hang, Valere Lambert, Eric Lindsey, Sagar Masuti, Takanori Matsuzawa, Jun Muto, Priyamvada Nanjundiah, Rino Salman, vSharadha Sathiakumar, and Harpreet Sethi.\njdpmoore/unicycle: Unicycle (Version 1.1). Zenodo.\nhttp://doi.org/10.5281/zenodo.5688288\n');
            if opts(1)
                fprintf('------\nBarbot S, Moore J D P, Lambert V.\nDisplacements and stress associated with distributed anelastic deformation in a half-space.\nBSSA (2017)\n');
            end
            if opts(2)
                fprintf('------\nMoore J D P and Sato D.\nThe role of topography, Earths curvature, and gravity in seismic cycles.\nIn prep (2022)\n');
            end
            if opts(3)
                fprintf('------\nLindsey et. al. 2019\nPlease email for reference (earth@jamesdpmoore.com)\n');
            end
            if opts(4)
                fprintf('------\nStress cycles\nPlease email for reference (earth@jamesdpmoore.com)\n');
            end
        end
        
    end
end

function [pos] = getPeriodicIndex(t,period)
% function GETPERIODICINDEX finds the positions in array t that
% are the closest to a periodic.
to=t(1);
pos=[];
for k=1:length(t)
    tc=t(k);
    if (tc-to>period)
        to=tc;
    end
    if tc>=to
        pos=[pos;k];
        to=tc+period;
    end
end
end
function [ix] = getClosestIndex(arr, value)
[~,ix] = min( abs( arr-value ) );
end