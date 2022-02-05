function [yp]= ode_acceleration_pwr(t,Y,o,style,stressing)
%
% Y = (s, tau, theta, V)
% yp = (V, dtau/dt, dtheta/dt, A)
%
yp=zeros(size(Y));

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       Fault                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
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

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Strain Volumes                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
if style == 1
    % deviatoric stress components
    s12p=Y(o.flt.N*o.flt.dgf+ 3:o.shz.dgf:end);
    s13p=Y(o.flt.N*o.flt.dgf+ 4:o.shz.dgf:end);
    
    %Kelvin self-stress tensor
    Q12k=s12p-o.shz.Gk.*Y(o.flt.N*o.flt.dgf+5:o.shz.dgf:end);
    Q13k=s13p-o.shz.Gk.*Y(o.flt.N*o.flt.dgf+6:o.shz.dgf:end);
    
    %Norm of the self-stress tensor
    QL2=sqrt(Q12k.^2+Q13k.^2);
    
    % Kelvin strain rates d epsilonK / dt = ( sigma' - 2G epsilonK ) / etaK
    e12k=Q12k.*(QL2.^(o.shz.m-1))./o.shz.etaK;
    e13k=Q13k.*(QL2.^(o.shz.m-1))./o.shz.etaK;
    
    % Norm of the shear stress tensor
    tauL2=sqrt(s12p.^2+s13p.^2);
    
    % Total strain = Kelvin + Maxwell strain rates
    e12=e12k+s12p.*(tauL2.^(o.shz.n-1))./o.shz.etaM;
    e13=e13k+s13p.*(tauL2.^(o.shz.n-1))./o.shz.etaM;
    
    % rate of strain
    yp(o.flt.N*o.flt.dgf+1:o.shz.dgf:end)=e12;
    yp(o.flt.N*o.flt.dgf+2:o.shz.dgf:end)=e13;
    
    % rate of Kelvin strain
    yp(o.flt.N*o.flt.dgf+5:o.shz.dgf:end)=e12k;
    yp(o.flt.N*o.flt.dgf+6:o.shz.dgf:end)=e13k;

end
if style == 2
    % isotropic stress
    p=(Y(o.flt.N*o.flt.dgf+4:o.shz.dgf:end)+Y(o.flt.N*o.flt.dgf+6:o.shz.dgf:end))/2;
    
    % deviatoric stress components
    s22p=Y(o.flt.N*o.flt.dgf+4:o.shz.dgf:end)-p;
    s23p=Y(o.flt.N*o.flt.dgf+5:o.shz.dgf:end);
    s33p=Y(o.flt.N*o.flt.dgf+6:o.shz.dgf:end)-p;
    
    % Maxwell strain rates
    e22=s22p./o.shz.etaM;
    e23=s23p./o.shz.etaM;
    e33=s33p./o.shz.etaM;
    
    % rate of strain
    yp(o.flt.N*o.flt.dgf+1:o.shz.dgf:end)=e22;
    yp(o.flt.N*o.flt.dgf+2:o.shz.dgf:end)=e23;
    yp(o.flt.N*o.flt.dgf+3:o.shz.dgf:end)=e33;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                  Stress Rates                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
if style == 1
    yp(2:o.flt.dgf:o.flt.N*o.flt.dgf) = acc - damping.*V ...
        +o.LK{2,style}*e12+o.LK{3,style}*e13;
    yp(o.flt.N*o.flt.dgf+ 3:o.shz.dgf:end)=o.LL{2,2}*e12+o.LL{3,2}*e13+o.KL{style,2}*(V-o.flt.Vpl);
    yp(o.flt.N*o.flt.dgf+ 4:o.shz.dgf:end)=o.LL{2,3}*e12+o.LL{3,3}*e13+o.KL{style,3}*(V-o.flt.Vpl);
end
if style == 2
    yp(2:o.flt.dgf:o.flt.N*o.flt.dgf) = acc - damping.*V ...
        +o.LK{4,style}*e22+o.LK{5,style}*e23+o.LK{6,style}*e33;
    yp(o.flt.N*o.flt.dgf+ 4:o.shz.dgf:end)=o.LL{4,4}*e22+o.LL{5,4}*e23+o.LL{6,4}*e33 ...
        +o.KL{style,4}*(V-o.flt.Vpl);
    yp(o.flt.N*o.flt.dgf+ 5:o.shz.dgf:end)=o.LL{4,5}*e22+o.LL{5,5}*e23+o.LL{6,5}*e33 ...
        +o.KL{style,5}*(V-o.flt.Vpl);
    yp(o.flt.N*o.flt.dgf+ 6:o.shz.dgf:end)=o.LL{4,6}*e22+o.LL{5,6}*e23+o.LL{6,6}*e33 ...
        +o.KL{style,6}*(V-o.flt.Vpl);
end
end