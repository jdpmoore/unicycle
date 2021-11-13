classdef coseismicPatch < unicycle.geometry.source
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
        % time to event
        t0;
        % name of event
        name;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = coseismicPatch(varargin)
            % PATCH is a class representing the geometry and physical
            % properties of fault patches, including position, orientation,
            % dimension and friction properties.
            %
            %   src = geometry.coseismicPatch('path/to/source/model.flt',t0, earthModel)
            %
            % creates a instance of fault patches with a slip distribution
            % defined in 'path/to/source/model.flt' to occur at time t0.
            % The stress components are computed using a specified Earth
            % model.
            % 
            % SEE ALSO: unicycle
            
            % project root directory
            obj@unicycle.geometry.source(varargin{1},varargin{3});
            
            if isempty(varargin)
                return
            end
            
            obj.name=varargin{1};
            obj.t0=varargin{2};
            
        end % constructor
        
        function tri=toCoseismicTriangle(obj)
            % TOCOSEISMICTRIANGLE converts the slip distribution to a
            % triangle mesh.
            %
            %   coseismicTriangle=obj.toCoseismicTriangle();
            %
            % INPUT:
            %   obj - instanced object of class coseismicPatch
            %
            % OUTPUT:
            %   coseismicTriangle - object of class coseismicTriangle
            %
            % SEE ALSO: unicycle
            
            import unicycle.geometry.triangle;
            
            tri=unicycle.geometry.coseismicTriangle();
            
            % time
            tri.t0=obj.t0;
            
            % number of triangles
            tri.N=2*obj.N;
            
            % identification number
            tri.id=1:tri.N;
            
            % triangular vertices
            tri.x=[ ...
                obj.x; ...
                obj.x+repmat(obj.L,1,3).*obj.sv; ...
                obj.x+repmat(obj.L,1,3).*obj.sv-repmat(obj.W,1,3).*obj.dv; ...
                obj.x-repmat(obj.W,1,3).*obj.dv];
            
            % mesh information (vertices of triangles)
            position=(1:obj.N)';
            tri.vertices=[ ...
                [position,obj.N+position,2*obj.N+position];
                [position,2*obj.N+position,3*obj.N+position]];
            
            % triangular center position
            tri.xc=[(tri.x(tri.vertices(:,1),1)+ ...
                     tri.x(tri.vertices(:,2),1)+ ...
                     tri.x(tri.vertices(:,3),1))/3, ...
                    (tri.x(tri.vertices(:,1),2)+ ...
                     tri.x(tri.vertices(:,2),2)+ ...
                     tri.x(tri.vertices(:,3),2))/3, ...
                    (tri.x(tri.vertices(:,1),3)+ ...
                     tri.x(tri.vertices(:,2),3)+ ...
                     tri.x(tri.vertices(:,3),3))/3];
            
            % unit vectors
            [tri.nv,tri.sv,tri.dv]=triangle.computeUnitVectors(tri.x,tri.vertices);
            
            % slip (amplitude, rake)
            tri.slip=[obj.slip;obj.slip];
            tri.rake=[obj.rake;obj.rake];
        end

    end % methods
   
end