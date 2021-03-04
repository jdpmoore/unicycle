 function [Yp]= odeViscoelastic_acceleration(~,Y,ss)            
% function odeViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%        | log(theta Vo / L) |           
%        |   log( V / Vo )   |
%    y = |       ...         |
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

% State variable
th   = Y(3:ss.dgfF:ss.M*ss.dgfF);

% Slip velocity 
V    = ss.Vo.*exp(Y(4:ss.dgfF:ss.M*ss.dgfF));

% Shear stress in zones of distributed deformation
tau12 = Y(ss.M*ss.dgfF+1:ss.dgfS:end);  
tau13 = Y(ss.M*ss.dgfF+2:ss.dgfS:end);
tau   = sqrt(tau12.^2+tau13.^2);

% Dislocation strain rate
Aeff = ss.Const_dis .* (tau).^(ss.n-1);
e12p = tau12 .* Aeff;
e13p = tau13 .* Aeff;

% Initiate state derivative
Yp=zeros(size(Y));  

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       Fault                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Slip velocity
Yp(1:ss.dgfF:ss.M*ss.dgfF) = V;

% Rate of state (rate of log(theta/theta0))
dth=(ss.Vo.*exp(-th)-V)./ss.L;
Yp(3:ss.dgfF:ss.M*ss.dgfF) = dth;

% Acceleration (rate of log(V/Vo))
kv = ss.K     *(V-ss.V_plate) + ...
     ss.k1212f*(e12p-ss.e12p_plate) + ...
     ss.k1312f*(e13p-ss.e13p_plate);
Yp(4:ss.dgfF:ss.M*ss.dgfF) = (kv - ss.b.*ss.sigmab.*dth)./ ...
                             (ss.a.*ss.sigmab + ss.damping.*V);

% Shear stress rate on fault due to fault and shear zones
Yp(2:ss.dgfF:ss.M*ss.dgfF) = kv - ss.damping.*V;

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     Shear Zones                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Stress rate due to shear zones and fault slip velocity
Yp(ss.M*ss.dgfF+1 : ss.dgfS : end) = ss.k1212*(e12p-ss.e12p_plate) + ...
                                     ss.k1312*(e13p-ss.e13p_plate) + ...
                                       ss.k12*(V-ss.V_plate);
                                   
Yp(ss.M*ss.dgfF+2 : ss.dgfS : end) = ss.k1213*(e12p-ss.e12p_plate) + ...
                                     ss.k1313*(e13p-ss.e13p_plate) + ...
                                       ss.k13*(V-ss.V_plate);

% Strain rate
Yp(ss.M*ss.dgfF+3 : ss.dgfS : end) = e12p;
Yp(ss.M*ss.dgfF+4 : ss.dgfS : end) = e13p;

end