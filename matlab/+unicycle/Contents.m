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
% Unicycle (Unified Cycle of Earthquakes) is a framework for the simulation
% of fault slip and distributed deformation using the integral method
% under the radiation damping approximation.
%
% Geometry
%   unicycle/geometry/coseismicPatch         - coseismic slip distribution
%   unicycle/geometry/patch                  - fault patch
%   unicycle/geometry/passiveReceiver        - passive stressed faults (no motion)
%   unicycle/geometry/source                 - moving fault loading receiver faults
%   unicycle/geometry/receiver               - receiver fault with slip evolution
%   unicycle/geometry/segment                - fault segment grouping patches
%   unicycle/geometry/triangle               - triangular fault patch
%   unicycle/geometry/triangleReceiver       - triangular receiver w/ slip evolution
%
% Green's function
%   unicycle/greens/stress                   - stress on receiver faults
%   unicycle/greens/stressKernels            - Okada (1992) stress kernels
%   unicycle/greens/triangleStressKernels    - Meade (2007) stress kernels
%   unicycle/greens/okada85                  - displacement at surface
%   unicycle/greens/shearZone16              - Distributed deformation using 
%                                              the analytical solution of James D P Moore
%                                              published in Barbot, Moore, and Lambert (2017)
%
%   hmmvp/unihmmvp                           - Green's functions with Hierarchical Matrices
%   
% Manifold
%   unicycle/manifold/gps                    - create a gps object
%
% Ordinary Differential Equation
%   unicycle/ode/evolution                   - models of fault slip evolution
%   unicycle/ode/rateandstate                - rate-and-state friction
%   unicycle/ode/rateandstatedamping         - rate-and-state friction w/ radiation damping
%   unicycle/ode/ratestrengthening           - rate-strengthening approximation
%   unicycle/ode/ratestrengthening_prestress - rate-strengthening friction w/ pre-stress
%
% Input/Output and formats
%   unicycle/export/exportflt_rfaults        - compatible w/ Relax, EDCMP, Gamra
%   unicycle/export/exportvtk_rfaults        - 3D visualization w/ Paraview
%   unicycle/export/exportxyz_rfaults        - GMT ASCII format
%   unicycle/export/grdread                  - read GMT GRD format
%
% Optimization
%   unicycle/optim/sim_anl                   - simulated annealing
%   unicycle/optim/mh                        - Metropolis Hastings
%   unicycle/optim/laplacian                 - smoothing matrix operator
%
% Authors:
% James D P Moore, Sylvain Barbot, 
% Lujia Feng, Yu Hang, Valere Lambert, Eric Lindsey, Sagar Masuti,       
% Takanori Matsuzawa, Jun Muto, Priyamvada Nanjundiah, Rino Salman,      
% Sharadha Sathiakumar, and Harpreet Sethi.
% Copyright 2017 Nanyang Technological University, Singapore
