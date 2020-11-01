% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                         %
%         E a r t h q u a k e   c y c l e s               %
%      o n   a   s t r i k e - s l i p   f a u l t        %
%                                                         %
% DESCRIPTION:                                            %
% Solves the governing equations for fault slip evolution % 
% using rate-and-state friction using the velocity in the %
% state vector and the acceleration in the state          %
% derivative.                                             %
%                                                         %
% Evaluates the evolution of slip using the integral      %
% method.                                                 %
%                                                         %
% AUTHOR:                                                 %
% Sylvain Barbot and Valere Lambert (April, 2017)         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all

minmax=@(x) [min(x(:)),max(x(:))];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            P H Y S I C A L   M O D E L                %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Rigidity (MPa)
G = 30e3;

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %

%                         0       
% _______________________________________________ 0 km   +--------> x2
%                         |                              |
%                         |                              |
%                         |      Brittle                 |
%                         |                              V
%                         |                              x3
%                         | fault        
% ----------------------------------------------- 35 km
%

% Fault Meshes
y3 = 0e3; % Fault starting depth (m)
y2 = 0e3; % Fault horizontal position (m)

% Brittle-Ductile Tranisition Depth (m)
Transition = 35e3;

% number of fault patches
ss.M = 120;

dz     = Transition/ss.M;
fpoles = y3+(0:ss.M)'*dz;

% top of fault patches
ss.y3f = fpoles(1:end-1); 
% width of fault patches
Wf     = ones(ss.M,1)*dz; 
  
%% Create stress kernels for fault interactions

ss.K=zeros(ss.M,ss.M);   % Fault self stress

% Evaluate the stress at the center of the fault
for k=1:ss.M
    % Stress on faults from fault slip
    ss.K(:,k)=s12h(y2,ss.y3f+dz/2,y2,ss.y3f(k),Wf(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% effective confining pressure on fault (MPa)
ss.sigmab = 1000;

% frictional parameters ( velocity-weakening friction, a-b < 0 )
ss.a = 1e-3*ones(size(ss.y3f));
ss.b = ss.a+2.5e-4*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.004*ones(size(ss.y3f));

% plate velocity (m/s)
ss.V_plate = 1e-9*ones(size(ss.y3f));

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(size(ss.y3f));

% shear wave speed (m/s)
ss.Vs = 3e3*ones(size(ss.y3f));

% radiation damping term
ss.damping=G./ss.Vs/2;

% Velocity-strengthening at top and bottom ( a-b > 0 )
% 5-15km velocity-weakening
top    = floor(0e3/(Transition/ss.M));
bottom = ceil(20e3/(Transition/ss.M));
ss.b(1:top)      = ss.a(1:top)     -2.1e-4*ones(top,1);
ss.b(bottom:end) = ss.a(bottom:end)-2.1e-4*ones(length(ss.a(bottom:end)),1);

% Fault Strength
ss.strength = ss.sigmab.*(ss.mu0+(ss.a-ss.b).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         N U M E R I C A L   S O L U T I O N          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% state parameters
ss.dgfF=4;

%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF,1);

% Fault patches
Y0(1:ss.dgfF:end) = zeros(size(ss.y3f));         % slip
Y0(2:ss.dgfF:end) = ss.strength;                 % stress
Y0(3:ss.dgfF:end) = log(ss.Vo./ss.V_plate);      % log(theta Vo / L)
Y0(4:ss.dgfF:end) = log(ss.V_plate*0.98./ss.Vo); % log(V/Vo)

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) earthquake_cycle_acceleration_powerlaw_ode(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6); 
[t,Y]=ode45(yp,[0 1e10],Y0,options);
toc

% Velocities
V=repmat(ss.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:ss.dgfF:ss.M*ss.dgfF));

% Maximum Velocity
Vmax=max(V,[],2);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

figure(2);clf;set(gcf,'name','Time Evolution')
f2a = subplot(3,1,1);cla;
pcolor(t(1:end-1)/3.15e7,ss.y3f/1e3,log10(V')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
colormap(f2a,parula);
title(h,'Slip Rate (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

f2b = subplot(3,1,2);cla;
plot(t(1:end-1)/3.15e7,log10(Vmax))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Maximum slip rate on fault')

f2c = subplot(3,1,3);cla;
plot(t(1:end-1)/3.15e7,log10(V(:,floor((top+bottom)/2))))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zone')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
figure(3);clf;set(gcf,'name','Time Step Evolution')

f3a=subplot(3,1,1);cla;
pcolor(1:length(t)-1,ss.y3f/1e3,log10(V')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
colormap(f3a,parula);
title(h,'Slip Rate (m/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

f3b = subplot(3,1,2);cla;
plot(1:length(t)-1,log10(Vmax))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Maximum slip rates on fault')

f3c = subplot(3,1,3);cla;
plot(1:length(t)-1,log10(V(:,floor((top+bottom)/2))))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zone')
