 function [Yp]= earthquake_cycle_acceleration_ode(~,Y,ss)            
% function odeViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is 
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        |   log( V / Vo )   |
%        \       ...         /
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

% Shear stress on faults
tauF = Y(2:ss.dgfF:end);

% State variables
th   = Y(3:ss.dgfF:end);

% Slip velocities 
V    = ss.Vo.*exp(Y(4:ss.dgfF:end));

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:ss.dgfF:ss.M*ss.dgfF) = V;

% Rate of state (rate of log(theta/theta0))
dth=(ss.Vo.*exp(-th)-V)./ss.L;
Yp(3:ss.dgfF:ss.M*ss.dgfF) = dth;

% Acceleration (rate of log(V/Vo))
kv = ss.K*(V - ss.V_plate);
Yp(4:ss.dgfF:ss.M*ss.dgfF) = (kv - ss.b.*ss.sigmab.*dth)./(ss.a.*ss.sigmab+ss.damping.*V);

% Shear stress rate on fault due to fault
Yp(2:ss.dgfF:ss.M*ss.dgfF) = kv - ss.damping.*V.*Yp(4:ss.dgfF:ss.M*ss.dgfF);

end

