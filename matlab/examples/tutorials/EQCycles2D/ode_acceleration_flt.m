function [yp]= ode_acceleration_flt(t,Y,o,style,stressing)
%
% Y = (s, tau, theta, V)
% Yp = (V, dtau/dt, dtheta/dt, A)
%
yp=zeros(size(Y));
theta = Y(3:o.flt.dgf:o.flt.N*o.flt.dgf); %state variable
V = o.flt.Vo.*exp(Y(4:o.flt.dgf:o.flt.N*o.flt.dgf)); % Slip velocity
if 0~=numel(o.flt.pinnedPosition) % pin patches
    theta(o.flt.pinnedPosition)=0;
    V(o.flt.pinnedPosition)=0;
end
yp(1:o.flt.dgf:o.flt.N*o.flt.dgf) = V; % Slip velocity
thetadot=(o.flt.Vo.*exp(-theta)-V)./o.flt.l; % Rate of state
yp(3:o.flt.dgf:o.flt.N*o.flt.dgf) = thetadot; % Rate of state
acc = o.KK{style,style}*(V-o.flt.Vpl); % Acceleration
sigmat=o.flt.sigma+stressing(t); % Normal stress
damping=0.5.*o.flt.earthModel.G./o.flt.Vs/2; %Radiation damping
yp(4:o.flt.dgf:o.flt.N*o.flt.dgf) = ((acc - o.flt.b.*sigmat.*thetadot)./ ...
    (o.flt.a.*sigmat + damping.*V)); % Acceleration
yp(2:o.flt.dgf:o.flt.N*o.flt.dgf) = acc - damping.*V; % Shear stress rate
end