classdef observationPoint < handle
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
        % index
        index;
        % position
        xc;
        % solution time
        t;
        % solution vector time series
        y;
        % number of time steps
        N;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = observationPoint(index)
            % OBSERVATIONPOINT is a class representing a patch and its
            % evolution.
            %
            %   opt = geometry.observationPoint(index)
            %
            % where index is the patch coordinate on a receiver fault.
            %
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle.
            
            obj.index=index;
            obj.N=0;
            
        end % constructor
        
    end % methods
    
end
