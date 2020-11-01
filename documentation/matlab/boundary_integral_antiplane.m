
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
% Evaluates the slip history on a fault in antiplane   %
% strain under the rate- and state-dependent friction  %
%                                                      %
% This code does not implement radiation damping and   %
% will fail to resolve quasi-dynamic ruptures. To      %
% resolve quasi-dynamic instabilities, use the code    %
% earthquake_cycle_acceleration.m                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            P H Y S I C A L   M O D E L               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% rigidity (MPa)
G=30e3;

% stress-interaction function
s12h=@(x2,x3,y2,y3,W) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-W)./((x2-y2).^2+(x3-y3-W).^2)-(x3+y3+W)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,W) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-W).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

M=30;
dz=0.15e3;
dx=0.3e3;
% top of slip patch
y3=(0:M-1)'*dz;
y2=(-M/2:M/2-1)*dx;
% width (down dip) of slip patch
W=ones(M,1)*dz;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%             S T R E S S   K E R N E L S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

ss.K=zeros(M,M);
for k=1:M
    % we evaluate the stress at the center 
    % of the slip patches
    ss.K(:,k)=s12h(0,y3+dz/2,0,y3(k),W(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% default friction properties (velocity-weakening friction)
% static friction coefficient
ss.mu0=0.1*ones(size(y3));
% frictional parameters
ss.a=1e-2*ones(size(y3));
ss.b=ss.a+1.4e-3*ones(size(y3));
% normal stress
ss.sigma=750.0*ones(size(y3));
% characteristic weakening distance
ss.L=0.1*ones(size(y3));
% plate velocity
ss.Vpl=1e-9*ones(size(y3));
% reference slip rate
ss.Vo=1e-6*ones(size(y3));

ss.b(1:4)=ss.a(1:4)-1.4e-3*ones(4,1);
ss.b(end-2:end)=ss.a(end-2:end)-1.4e-3*ones(3,1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         N U M E R I C A L   S O L U T I O N          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

ss.dgf=3;
Y0=zeros(M*ss.dgf,1);
Y0(1:ss.dgf:end)=ss.Vpl;
Y0(2:ss.dgf:end)=ss.mu0.*ss.sigma;
Y0(3:ss.dgf:end)=zeros(M,1);

% initialize the function handle with 
% set constitutive parameters
yp=@(t,y) odefunAntiplane(t,y,ss);

% solve the system
options=odeset('Refine',1,'RelTol',3e-14,'InitialStep',1e-7);
[t,Y]=ode23(yp,[0 3e7],Y0,options);

% compute the instantaneous derivative
Yp=zeros(size(Y));
for i=1:length(t)
    Yp(i,:)=yp(t(i),Y(i,:)');
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

figure(1);clf;set(gcf,'name','Time evolution')

subplot(2,1,1);cla;
pcolor(t/3.15e7,y3/1e3,log10(Yp(:,1:ss.dgf:end)')), shading flat
h=colorbar('Location','NorthOutside');
title(h,'Velocity (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

subplot(2,1,2);cla;
plot(t/3.15e7,log10(Yp(:,(M/2-1)*ss.dgf+1)'))
xlabel('Time (yr)')
ylabel('Velocity (m/s) log10')
title('Time series at the fault center')


figure(2);clf;set(gcf,'name','Evolution with time steps')

subplot(2,1,1);cla;
pcolor((1:length(t)),y3/1e3,log10(Yp(:,1:ss.dgf:end)')), shading flat
h=colorbar('Location','NorthOutside');
title(h,'Velocity (m/s)')
xlabel('Step')
ylabel('Depth (km)');

subplot(2,1,2);cla;
plot((1:length(t)),log10(Yp(:,(M/2-1)*ss.dgf+1)'))
xlabel('Step')
ylabel('Velocity (m/s) log10')
title('Evolution at the fault center')



