classdef gpsReceiver < unicycle.manifold.gps
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
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = gpsReceiver(network,evl,vecsize,obsname)
            
            
            % GPS is a class representing GPS data and Green's functions.
            %
            %   gps = manifold.gpsReceiver(network,evl,vecsize)
            %
            % creates a instance of GPS data connected to a simulation.
            %
            % INPUT:
            %
            % network
            %
            % if network is a filename, for example 'sopac.dat', containing
            %
            %          # i NAME x1 x2 x3
            %            1 GPS1  0  0  0
            %            2 GPS2  1  1  0
            %
            % alternatively network may be a vector of size N x vecsize
            % containing locations for evaluating displacements
            %
            % evl      object of type ode.evolution containing a time
            %          series of fault displacement
            % vecsize  number of components of displacement vector
            %
            % SEE ALSO: unicycle
            
            import unicycle.greens.*
            
            obj.vecsize=vecsize;
            
            if isnumeric(network)
                obj.x=network;
                obj.D=size(network,1);
                obj.stationName=1:obj.D;
            else
                filename=network;
                [~,obj.stationName,x1,x2,x3]=...
                    textread(filename,'%d %s %f %f %f','commentstyle','shell');
                obj.x=[x2,x1,0*x3];
                obj.D=size(x2,1);
            end
            
            
            %             if ~isempty(varargin)
            %                 if isempty(varargin{:})
            %                     knlpathobs=evl.knlpath;
            %                 else
            %                     % kernel path
            %                     o.knlpath=char(varargin{:});
            %                     if ~exist(o.knlpath,'dir')
            %                         st=strsplit(o.knlpath,'/');
            %                         mkdir(st{end-1});
            %                     end
            %                 end
            %             end
            % source Green's functions (includes strike slip and dip slip)
            
            if isobject(evl.src)
                obj.knl.FO={'s', 'd'};
                
                if ~exist(evl.knlpath)
                    [obj.FO{1,1},obj.FO{2,1}]=evl.src.displacementKernels(obj.x,vecsize);
                else
                    if ~exist(strcat(evl.knlpath,obsname,'FO_s.grd'),'file')
                        if ~exist(strcat(evl.knlpath,obsname,'FO.mat'),'file')
                            [obj.FO{1,1},obj.FO{2,1}]=evl.src.displacementKernels(obj.x,vecsize);
                            for i = 1:numel(obj.FO)
                                fname=strcat(evl.knlpath,obsname,'FO_',obj.knl.FO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.FO{i},fname)
                            end
                        else
                            load(strcat(evl.knlpath,obsname,'FO.mat'));
                            obj.FO=sFO;
                            clear sFO
                            for i = 1:numel(obj.FO)
                                fname=strcat(evl.knlpath,obsname,'FO_',obj.knl.FO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.FO{i},fname)
                            end
                        end
                    else
                        for i = 1:numel(obj.FO)
                            fname=strcat(evl.knlpath,obsname,'FO_',obj.knl.FO{i},'.grd');
                            [~,~,obj.FO{i}]=unicycle.export.grdread(fname);
                        end
                    end
                end
            end
            
            % receiver Green's functions (includes strike slip and dip slip)
            if isobject(evl.flt)
                obj.knl.KO={'s', 'd'};
                if ~exist(evl.knlpath)
                    [obj.KO{1,1},obj.KO{2,1}]=evl.flt.displacementKernels(obj.x,vecsize);
                else
                    if ~exist(strcat(evl.knlpath,obsname,'KO_s.grd'),'file')
                        if ~exist(strcat(evl.knlpath,obsname,'KO.mat'),'file')
                            [obj.KO{1,1},obj.KO{2,1}]=evl.flt.displacementKernels(obj.x,vecsize);
                            for i = 1:numel(obj.KO)
                                fname=strcat(evl.knlpath,obsname,'KO_',obj.knl.KO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.KO{i},fname)
                            end
                        else
                            load(strcat(evl.knlpath,obsname,'KO.mat'));
                            obj.KO=sKO;
                            clear sKO
                            for i = 1:numel(obj.KO)
                                fname=strcat(evl.knlpath,obsname,'KO_',obj.knl.KO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.KO{i},fname)
                            end
                        end
                    else
                        for i = 1:numel(obj.KO)
                            fname=strcat(evl.knlpath,obsname,'KO_',obj.knl.KO{i},'.grd');
                            [~,~,obj.KO{i}]=unicycle.export.grdread(fname);
                        end
                    end
                end
            end
            
            % Strain Volume Green's functions (includes strike slip and dip slip)
            if isobject(evl.shz)
                obj.knl.LO={'11','12','13','22','23','33'};
                if ~exist(evl.knlpath)
                    [obj.LO{1,1},obj.LO{2,1},obj.LO{3,1}, ...
                        obj.LO{4,1},obj.LO{5,1},obj.LO{6,1}]=evl.shz.displacementKernels(obj.x,vecsize);
                else
                    if ~exist(strcat(evl.knlpath,obsname,'LO_11.grd'),'file')
                        if ~exist(strcat(evl.knlpath,obsname,'LO.mat'),'file')
                            [obj.LO{1,1},obj.LO{2,1},obj.LO{3,1}, ...
                                obj.LO{4,1},obj.LO{5,1},obj.LO{6,1}]=evl.shz.displacementKernels(obj.x,vecsize);
                            for i = 1:numel(obj.LO)
                                fname=strcat(evl.knlpath,obsname,'LO_',obj.knl.LO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.LO{i},fname)
                            end
                        else
                            load(evl.knlpath,obsname,'LO.mat');
                            obj.LO=sLO;
                            clear sLO
                            for i = 1:numel(obj.LO)
                                fname=strcat(evl.knlpath,obsname,'LO_',obj.knl.LO{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],obj.LO{i},fname)
                            end
                        end
                    else
                        for i = 1:numel(obj.LO)
                            fname=strcat(evl.knlpath,obsname,'LO_',obj.knl.LO{i},'.grd');
                            [~,~,obj.LO{i}]=unicycle.export.grdread(fname);
                        end
                    end
                end
            end
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                obj.simulation(evl);
            end
            
        end % constructor
        
        
    end % methods
    
end
