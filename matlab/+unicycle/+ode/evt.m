classdef evt < handle
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
    % http://doi.org/10.5281/zenodo.5688288                                  %
    %                                                                        %
    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    properties
        % start and end times
        tStart;
        tEnd;
        
        % start and end indices
        iStart;
        iEnd;
        
        % solution vector at start and end times
        yStart;
        yEnd;
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = evt(tStart,tEnd,iStart,iEnd,yStart,yEnd)
            % EVENT is a class representing a slip event
            %
            %   opt = ode.event(tStart,tEnd,iStart,iEnd,yStart,yEnd)
            %
            % where xStart and xEnd and the start and end x attributes.
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle.
            
            obj.tStart=tStart;
            obj.tEnd=tEnd;
            
            obj.iStart=iStart;
            obj.iEnd=iEnd;
            
            obj.yStart=yStart;
            obj.yEnd=yEnd;
            
        end % constructor
        
    end % methods
    
end
