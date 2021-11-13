classdef gps < unicycle.manifold.gps
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
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = gps(filename,evl,vecsize,prefix)
            % GPS is a class representing GPS data and Green's functions.
            %
            %   gps = unicycle.manifold.edcmp.gps(network,evl,vecsize,prefix)
            %
            % creates a instance of GPS data.
            %
            % INPUT:
            %
            % network  filename, for example 'sopac.dat', containing
            %
            %          # i NAME x1 x2 x3
            %            1 GPS1  0  0  0
            %            2 GPS2  1  1  0
            %
            % evl      object of type ode.evolution containing a time
            %          series of fault displacement
            % vecsize  number of components of displacement vector
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            import unicycle.greens.edcmp
            
            obj.vecsize=vecsize;
            
            [~,obj.stationName,x1,x2,x3]=...
                textread(filename,'%d %s %f %f %f','commentstyle','shell');
            obj.x=[x2,x1,0*x3];
            obj.D=size(x2,1);
            
            % source Green's functions (includes strike slip and dip slip)
            obj.G=edcmp.G(evl.src,prefix,obj.x,vecsize);
            % receiver Green's functions (includes strike slip and dip slip)
            obj.H=edcmp.G(evl.rcv,prefix,obj.x,vecsize);
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                obj.simulation(evl);
            end
            
        end % constructor
        
    end % methods
    
end
