function [yp] = odefun(~,y,ss)
% function ODEFUN describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state 
% vector y is 
%
%        /        s          \
%    y = |       tau         |
%        \ log(theta Vo / L) /
%
% Instead of integrating numerically the againg law
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

V=2*ss.Vo*exp(...
    (y(2)-ss.mu0*ss.sigma-ss.b*ss.sigma*y(3))/(ss.a*ss.sigma));

yp=[V; ...
    ss.k*(ss.Vpl-V); ...
    (ss.Vo*exp(-y(3))-V)/ss.L];

end
