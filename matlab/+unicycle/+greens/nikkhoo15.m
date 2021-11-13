classdef nikkhoo15 < unicycle.greens.earthModel
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
        % rigidity
        G;
        % Poisson's ratio
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=nikkhoo15(G,nu)
            % NIKKHOO15 is a class providing the necessary functions to
            % compute the stress interactions between sources and receivers
            % and among the receivers using triangular dislocations.
            %
            %   earthModel = greens.nikkhoo15(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space.
            %
            % SEE ALSO: unicycle, geometry.triangle, geometry.triangleReceiver
            
            if (0==nargin)
                return
            end
            
            assert(0<=G,'unicycle.greens.nikkhoo15::rigidity must be positive.')
            assert(nu<=0.5,'unicycle.greens.nikkhoo15::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'unicycle.greens.nikkhoo15::Poisson''s ratio should be greater than -1.')
            
            o.G=G;
            o.nu=nu;
        end
        
        function [varargout]=tractionKernels(obj,src,rcv)
            % TRACTIONKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % rcv - receiver fault
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeTractionKernelsNikkhoo15(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=stressKernels(obj,src,rcv)
            % STRESSKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % rcv - point
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeStressKernelsNikkhoo15(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=displacementKernels(obj,src,x,vecsize)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % src - source fault
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeDisplacementKernelsNikkhoo15(src,obj.nu,x,vecsize);
        end
    end % methods
end % class definition
