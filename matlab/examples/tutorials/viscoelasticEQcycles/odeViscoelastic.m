 function [Yp]= odeViscoelastic(~,Y,ss)            
% function odeViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%    y = |       tau         |         
%        | log(theta Vo / L) |           
%        |       ...         |
%        |       s12         |
%        |       s13         |
%        |       e12         |
%        \       e13         /
%
% Instead of integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% Considering the following relationships
%
%    d phi = d theta / theta
%
%    d theta = d phi theta = d phi exp(phi) L / Vo
%
%    d theta / dt = d phi exp(phi) / dt L / Vo
%
%    d phi exp(phi) / dt L / Vo = 1 - V theta / L = 1 - V exp(phi) / Vo
%
% we obtain the evolution law for the new variable
% 
%    d phi / dt = ( Vo exp(-phi) - V ) / L
%
%  The state vector is split with initial cells considering fault parameters
%  and the latter portion considering strain in finite volumes with (s12,s13)
%  and (e12,e13) being the 2D antiplane shear stress and strain components.

G = 30e3; % MPa

% Shear stress on faults
tauF_W = Y(2:ss.dgfF:ss.M*ss.dgfF);
tauF_E = Y(ss.M*ss.dgfF+2:ss.dgfF:2*ss.M*ss.dgfF);

% State variables
th_W   = Y(3:ss.dgfF:ss.M*ss.dgfF);
th_E   = Y(ss.M*ss.dgfF+3:ss.dgfF:2*ss.M*ss.dgfF);

% Slip velocities 
% Use the Lambert W function to approximate the slip velocity incorporating
% radiation damping
V_W = (2.*ss.Vs.*ss.aW.*ss.sigmab./G).*...
      Lambert_W(G*ss.Vo./(2*ss.Vs.*ss.aW.*ss.sigmab).*...
      exp((tauF_W-ss.mu0.*ss.sigmab-ss.sigmab.*ss.bW.*th_W)./(ss.sigmab.*ss.aW)));

V_E = (2.*ss.Vs.*ss.aE.*ss.sigmab./G).*...
      Lambert_W(G*ss.Vo./(2*ss.Vs.*ss.aE.*ss.sigmab).*...
      exp((tauF_E-ss.mu0.*ss.sigmab-ss.sigmab.*ss.bE.*th_E)./(ss.sigmab.*ss.aE)));

% Shear stress in zones of distributed deformation
tau12 = Y(2*ss.M*ss.dgfF+1:ss.dgfS:end);  
tau13 = Y(2*ss.M*ss.dgfF+2:ss.dgfS:end);
tau   = sqrt(tau12.^2+tau13.^2);

% Dislocation strain rate
Aeff = ss.Const_dis .* (tau).^(ss.n-1);
e12p = tau12 .* Aeff;
e13p = tau13 .* Aeff;

% Initiate state derivative
Yp=zeros(size(Y));  

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     West Fault                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Slip velocity
Yp(1:ss.dgfF:ss.M*ss.dgfF) = V_W;

% Shear stress rate on fault due to fault and shear zones
Yp(2:ss.dgfF:ss.M*ss.dgfF) = ss.KWW    *(V_W-ss.V_plate)     + ss.KEW    *(V_E-ss.V_plate)+...
                             ss.k1212fW*(e12p-ss.e12p_plate) + ss.k1312fW*(e13p-ss.e13p_plate);

% Rate of state
Yp(3:ss.dgfF:ss.M*ss.dgfF) = (ss.Vo.*exp(-th_W)-V_W)./ss.L;

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     East Fault                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Slip velocity
Yp(ss.M*ss.dgfF+1:ss.dgfF:2*ss.M*ss.dgfF) = V_E;

% Shear stress rate on fault due to fault and shear zones
Yp(ss.M*ss.dgfF+2:ss.dgfF:2*ss.M*ss.dgfF) = ss.KEE    *(V_E-ss.V_plate)     + ss.KWE    *(V_W-ss.V_plate) +...
                                            ss.k1212fE*(e12p-ss.e12p_plate) + ss.k1312fE*(e13p-ss.e13p_plate);

% Rate of state
Yp(ss.M*ss.dgfF+3:ss.dgfF:2*ss.M*ss.dgfF) = (ss.Vo.*exp(-th_E)-V_E)./ss.L;

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     Shear Zones                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Stress rate due to shear zones and fault slip velocity
Yp(2*ss.M*ss.dgfF+1 : ss.dgfS : end) = ss.k1212*(e12p-ss.e12p_plate) + ss.k1312*(e13p-ss.e13p_plate) + ...
                                       ss.k12E *(V_E-ss.V_plate)     + ss.k12W *(V_W-ss.V_plate);
                                   
Yp(2*ss.M*ss.dgfF+2 : ss.dgfS : end) = ss.k1213*(e12p-ss.e12p_plate) + ss.k1313*(e13p-ss.e13p_plate) + ...
                                       ss.k13E *(V_E-ss.V_plate)     + ss.k13W *(V_W-ss.V_plate);

% Strain rate
Yp(2*ss.M*ss.dgfF+3 : ss.dgfS : end) = e12p;
Yp(2*ss.M*ss.dgfF+4 : ss.dgfS : end) = e13p;


end