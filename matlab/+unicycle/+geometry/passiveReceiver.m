classdef passiveReceiver < unicycle.geometry.patch
    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Copyright 2017 Nanyang Technological University, Singapore             %
    %                                                                        %
    % This file is part of UNICYCLE                                          %
    %                                                                        %
    % UNICYCLE is free software for non commercial usage:                    %
    % you can redistribute it and/or modify it under the terms of the        %
    % Creative Commons CC BY-NC-SA 4.0 licence, please see                   %
    % https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode            %
    %                                                                        %
    % UNICYCLE is distributed in the hope that it will be useful,            %
    % but WITHOUT ANY WARRANTY; without even the implied warranty of         %
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   %
    %                                                                        %
    % All intellectual property rights and licences for commercial usage     %
    % are reserved by Nanyang Technological University, please contact       %
    % NTUItive if you are interested in a commericial licence for UNICYCLE.  %
    % https://www.ntuitive.sg/                                               %
    %                                                                        %
    % If you use this code, please cite as James D P Moore, Sylvain Barbot,  %
    % Lujia Feng, Yu Hang, Valere Lambert, Eric Lindsey, Sagar Masuti,       %
    % Takanori Matsuzawa, Jun Muto, Priyamvada Nanjundiah, Rino Salman,      %
    % Sharadha Sathiakumar, and Harpreet Sethi. (2019, September 25).        %
    % jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo.                     %
    % http://doi.org/10.5281/zenodo.4471162                                  %
    %                                                                        %
    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    properties
        % static friction
        mu0;
        
		% shear stress
		% KK{i,j} shear stress s{i} due to fault slip d{i} where
		% s={strike shear, dip shear, normal stress} and d={strike slip, dip slip}
		KK=cell(2,3);
		
        % modeled shear stress in strike direction due to rcv and source
        sr,ss;
        % modeled shear stress in dip direction due to receiver and source
        dr,ds;
        
        % modeled tension due to receiver and source
        tr,ts;
        % modeled pressure due to receiver and source
        pr,ps;
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = passiveReceiver(basename,evl)
            % PASSIVERECEIVER is a class representing the geometry and 
            % physicalproperties of receiver fault patches, including 
            % position, orientation, dimension and friction properties.
            %
            %   src = geometry.passivereceiver('basename')
            %
            % creates a instance of fault patches for a passive receiver 
            % fault for stress change calculations.
            %
            % SEE ALSO: unicycle
            
            import unicycle.greens.*
            
            fname=[basename '_seg.flt'];
            if 2==exist(fname,'file')
                [fm,seg]=obj.seg2flt(fname);
            else
                fname=[basename '_patch.flt'];
                assert(2==exist(fname,'file'));
                [id,x1,x2,x3,len,width,str,d,rak]=...
                    textread(fname,'%u %f %f %f %f %f %f %f %f',...
                    'commentstyle','shell');
                obj.id=id;
                fm=[x1,x2,x3,len,width,str,d,rak];
                seg={unicycle.geometry.segment(mean(x2),mean(x1),mean(-x3),...
                    mean(len),mean(width),mean(str),mean(d),mean(rak),...
                    0,length(x1))};
            end
            
            % patch properties
            obj.N=size(fm,1);
            obj.slip=zeros(obj.N,1);
            obj.x=[fm(:,[2,1]),-fm(:,3)];
            obj.L=fm(:,4);
            obj.W=fm(:,5);
            obj.strike=fm(:,6);
            obj.dip=fm(:,7);
            obj.rake=fm(:,8);
            
            % default static friction
            obj.mu0=obj.L*0+0.6;
            
            % segments
            obj.segments=seg;
            
            % unit vectors in the strike direction
            obj.sv=[sind(obj.strike),cosd(obj.strike),zeros(obj.N,1)];
            
            % unit vectors in the dip direction
            obj.dv=[...
                -cosd(obj.strike).*cosd(obj.dip), ...
                +sind(obj.strike).*cosd(obj.dip), ...
                +sind(obj.dip)];
            
            % unit vectors in the normal direction
            obj.nv=[...
                +cosd(obj.strike).*sind(obj.dip), ...
                -sind(obj.strike).*sind(obj.dip), ...
                +cosd(obj.dip)];
            
            % center of fault patch
            obj.xc=[...
                obj.x(:,1)+obj.L/2.*obj.sv(:,1)-obj.W/2.*obj.dv(:,1),...
                obj.x(:,2)+obj.L/2.*obj.sv(:,2)-obj.W/2.*obj.dv(:,2),...
                obj.x(:,3)+obj.L/2.*obj.sv(:,3)-obj.W/2.*obj.dv(:,3)];
            
            % receiver stress interactions including pressure
            %[obj.Kss,obj.Ksd,obj.Kds,obj.Kdd, obj.SS,obj.DS]=evl.stressKernels(obj);
            [obj.KK{:}]=evl.flt.stressKernels(obj);
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                obj.sr=zeros(obj.N,numel(evl.t));
                obj.dr=zeros(obj.N,numel(evl.t));
                obj.tr=zeros(obj.N,numel(evl.t));
                obj.pr=zeros(obj.N,numel(evl.t));
                
                for k=1:length(evl.t)
                    % receiver shear stress in strike direction model
                    %obj.sr(:,k)=obj.Kss*evl.y(1:evl.flt.dgf:end,k)+...
                    %            obj.Kds*evl.y(2:evl.flt.dgf:end,k);
                   
                   
                    obj.sr(:,k)=obj.KK{1}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                                obj.KK{3}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                    % receiver pressure model due to strike slip and dip slip
                    %obj.dr(:,k)=obj.Ksd*evl.y(1:evl.flt.dgf:end,k)+...
                    %            obj.Kdd*evl.y(2:evl.flt.dgf:end,k);
                    obj.dr(:,k)=obj.KK{2}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                                obj.KK{4}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                    
                    % receiver tension model due to strike slip and dip slip
                    %obj.tr(:,k)=obj.SS.tn*evl.y(1:evl.flt.dgf:end,k)+...
                    %            obj.DS.tn*evl.y(2:evl.flt.dgf:end,k);
                    obj.tr(:,k)=obj.KK{5}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                                obj.KK{6}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                    % receiver pressure model due to strike slip and dip slip
                    %obj.pr(:,k)=-(obj.SS.Kxx+obj.SS.Kyy+obj.SS.Kzz)*evl.y(1:evl.flt.dgf:end,k) ...
                    %           -(obj.DS.Kxx+obj.DS.Kyy+obj.DS.Kzz)*evl.y(2:evl.flt.dgf:end,k);
                    obj.pr(:,k)=-(obj.KK{5})*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k) ...
                                -(obj.KK{6})*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                end
            end
        end % constructor
        
        
        function [] = simulation(obj,evl)
            % method SIMULATION builds forward models of stress if
            % simulation exists
            %
            % SEE ALSO: unicycle
            
            obj.sr=zeros(obj.N,numel(evl.t));
            obj.dr=zeros(obj.N,numel(evl.t));
            obj.tr=zeros(obj.N,numel(evl.t));
            obj.pr=zeros(obj.N,numel(evl.t));
            
            for k=1:length(evl.t)
                
                % receiver shear stress in strike direction model
                obj.sr(:,k)=obj.KK{1}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                            obj.KK{3}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                % receiver pressure model due to strike slip and dip slip
                obj.dr(:,k)=obj.KK{2}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                            obj.KK{4}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                        
                % receiver tension model due to strike slip and dip slip
                obj.tr(:,k)=obj.KK{5}*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                            obj.KK{6}*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
                % receiver pressure model due to strike slip and dip slip
                obj.pr(:,k)=-(obj.KK{5})*evl.y(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k)+...
                            -(obk.KK{6})*evl.y(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N,k);
            end
        end % simulation
        
        function [tpr] = toTrianglePassiveReceiver(obj,evl)
            % TOTRIANGLEPASSIVERECEIVER provides a triangle representation
            % of the passive receiver.
            %
            %   obj.toTrianglePassiveReceiver(evl)
            %
            % SEE ALSO: unicycle
            
            tpr=unicycle.geometry.trianglePassiveReceiver();
            tpr=obj.toTriangle(tpr);
            
            % receiver stress interactions including pressure
            [tpr.Kss,tpr.Ksd,tpr.Kds,tpr.Kdd, ...
                tpr.SS,tpr.DS]=evl.flt.stressKernels(tpr);
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                tpr.sr=zeros(tpr.N,numel(evl.t));
                tpr.dr=zeros(tpr.N,numel(evl.t));
                tpr.tr=zeros(tpr.N,numel(evl.t));
                tpr.pr=zeros(tpr.N,numel(evl.t));
                
                for k=1:length(evl.t)
                    % receiver shear stress in strike direction model
                    tpr.sr(:,k)=tpr.Kss*evl.y(1:evl.dgf:end,k)+...
                                tpr.Kds*evl.y(2:evl.dgf:end,k);
                    % receiver pressure model due to strike slip and dip slip
                    tpr.dr(:,k)=tpr.Ksd*evl.y(1:evl.dgf:end,k)+...
                                tpr.Kdd*evl.y(2:evl.dgf:end,k);

                    % receiver tension model due to strike slip and dip slip
                    tpr.tr(:,k)=tpr.SS.Ktn*evl.y(1:evl.dgf:end,k)+...
                                tpr.DS.Ktn*evl.y(2:evl.dgf:end,k);
                    % receiver pressure model due to strike slip and dip slip
                    tpr.pr(:,k)=-(tpr.SS.Kxx+tpr.SS.Kyy+tpr.SS.Kzz)*evl.y(1:evl.dgf:end,k) ...
                                -(tpr.DS.Kxx+tpr.DS.Kyy+tpr.DS.Kzz)*evl.y(2:evl.dgf:end,k);
                end
            end
        end
        
    end % methods
    
end
