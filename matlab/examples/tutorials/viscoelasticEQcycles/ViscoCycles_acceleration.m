% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            Unified earthquake cycles of               %
%          fault slip and viscoelastic strain           %
%                                                       %
% DESCRIPTION:                                          %
% Simulates unified earthquake cycles including slip on %
% on a strike-slip fault and viscoelastic strain below  %
% the brittle-ductile transition in 2D antiplane.       %
%                                                       %
% Solves the governing equations for fault slip         %
% evolution using rate-and-state friction coupled       %
% with distributed viscoelastic deformation governed by %
% dislocation creep.                                    %
%                                                       %
% Utilizes solutions for displacement and stress fields %
% due to slip and distributed strain in 2D to couple    %
% fault and off-fault deformation. Evaluates the        %
% evolution of slip and strain using the integral       %
% method.                                               %
%                                                       %
% AUTHORS:                                               %
% Valere Lambert and Sylvain Barbot (April, 2017)       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all

% Minimum and maximum function
minmax=@(x) [min(x(:)),max(x(:))];

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            P H Y S I C A L   M O D E L                %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Rigidity (MPa)
G = 30e3;

% Stress kernels for distributed deformation
s1312=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    +log((x2-L/2).^2+(x3+D+W).^2) - log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    -log((x2-L/2).^2+(x3+D).^2) + log((x2+L/2).^2+(x3+D).^2));

s1212=@(D,L,W,x2,x3) G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2)) ...
    -atan((x3+D+W)./(x2-L/2))-atan((x3+D)./(x2+L/2)) ...
    +atan((x3+D)./(x2-L/2))+atan((x3+D+W)./(x2+L/2)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

s1213=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3+D+W).^2) + log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    +log((x2-L/2).^2+(x3+D).^2) - log((x2+L/2).^2+(x3+D).^2));

s1313=@(D,L,W,x2,x3) G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W)) ...
    +atan((x2+L/2)./(x3+D))  -atan((x2-L/2)./(x3+D)) ...
    -atan((x2+L/2)./(x3+D+W))+atan((x2-L/2)./(x3+D+W)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Displacement Kernels for distributed deformation
uk12=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x3-D-W).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x3-D-W).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x3-D).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x3-D).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x2-L/2).*(atan((x3-D-W)./(x2-L/2))-atan((x3-D)./(x2-L/2))) ...
    +2*(x2+L/2).*(atan((x3-D)./(x2+L/2))-atan((x3-D-W)./(x2+L/2))) ...
    +(x3+D+W).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x3+D+W).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x3+D).*log((x2+L/2).^2+(x3+D).^2) ...
    +(x3+D).*log((x2-L/2).^2+(x3+D).^2) ...
    +2*(x2+L/2).*(atan((x3+D+W)./(x2+L/2))-atan((x3+D)./(x2+L/2))) ...
    +2*(x2-L/2).*(atan((x3+D)./(x2-L/2))-atan((x3+D+W)./(x2-L/2))) ...
    );

uk13=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x2-L/2).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x3-W-D).*(atan((x2-L/2)./(x3-D-W))-atan((x2+L/2)./(x3-D-W))) ...
    +2*(x3-D).*(atan((x2+L/2)./(x3-D))-atan((x2-L/2)./(x3-D))) ...
    +(x2-L/2).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3+D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3+D).^2) ...
    +2*(x3+W+D).*(atan((x2-L/2)./(x3+D+W))-atan((x2+L/2)./(x3+D+W))) ...
    +2*(x3+D).*(atan((x2+L/2)./(x3+D))-atan((x2-L/2)./(x3+D))) ...
    );

% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
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
%              Viscoelastic Shear Zones
%

% Fault Meshes
y3 = 0e3; % Faul starting depth (m)
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

%% Shear Zone Mesh
ss.Nx=50;
ss.Nz=50;
eps=1e-12;

% Patch edges along x3 (Width increases with depth)
ss.polesz = Transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*Transition;
% center of shear zone (x3)
ss.x3c    = (ss.polesz(2:end)+ss.polesz(1:end-1))/2; 
% shear zone width 
W         = ss.polesz(2:end)-ss.polesz(1:end-1);     

xx3t = repmat(ss.polesz(1:end-1)',ss.Nx,1);  % tops
xx3b = repmat(ss.polesz(2:end)',  ss.Nx,1);  % bottoms
xx3c = repmat(ss.x3c',ss.Nx,1);              % centers

% Patch edges along x2
% Make boxes of uniform length (L) from -50km to 50km
polesxc = (-50:2:50)'*1e3;

% Length increases with distance from x2=0 outside of window
edges   = floor(ss.Nx-length(polesxc)+1)/2; 
polesxl = flipud(min(polesxc)-tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition);
polesxr = max(polesxc)+tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition;

ss.polesx = [polesxl(1:end-1);polesxc;polesxr(2:end)];
% center of shear zone (x2)
ss.x2c    = (ss.polesx(2:end)+ss.polesx(1:end-1))/2;     
% shear zone length 
L         = ss.polesx(2:end)-ss.polesx(1:end-1);             

xx2l = repmat(ss.polesx(1:end-1),1,ss.Nz);  % left edge
xx2r = repmat(ss.polesx(2:end),  1,ss.Nz);  % right edge
xx2c = repmat(ss.x2c,1,ss.Nz);              % center

%% Grid points for plotting patches
xpointsF = [repmat(y2-0.15e3,length(ss.y3f),1),...
            repmat(y2+0.15e3,length(ss.y3f),1),...
            repmat(y2+0.15e3,length(ss.y3f),1),...
            repmat(y2-0.15e3,length(ss.y3f),1)]';
      
zpointsF = [repmat(fpoles(1:end-1),1,1),...
            repmat(fpoles(1:end-1),1,1),...
            repmat(fpoles(2:end),1,1),...
            repmat(fpoles(2:end),1,1)]';

xpoints = [xx2l(:)';xx2r(:)';xx2r(:)';xx2l(:)'];
zpoints = [xx3t(:)';xx3t(:)';xx3b(:)';xx3b(:)'];
  
%% GPS points for surface displacement
x2GPS = (-400:10:399)'*1e3;

%% Create stress kernels for fault and shear zone interactions

% Stress kernels from fault
ss.k12=zeros(length(xx2c(:)),ss.M);
ss.k13=zeros(length(xx2c(:)),ss.M);

ss.K=zeros(ss.M,ss.M);   % Fault self stress

% Displacement kernels
ss.ku1=zeros(length(x2GPS),ss.M);

% Fields from Faults
% Evaluate the stress at the center of the fault and shear zone patches
% and displacement at the top of the patch
for k=1:ss.M
    % Stress on shear zones from fault slip
    ss.k12(:,k)=s12h(xx2c(:),xx3c(:),y2,ss.y3f(k),Wf(k));  
    ss.k13(:,k)=s13h(xx2c(:),xx3c(:),y2,ss.y3f(k),Wf(k));
        
    % Stress on faults from fault slip
    ss.K(:,k)=s12h(y2,ss.y3f+dz/2,y2,ss.y3f(k),Wf(k));
    
    % Dispacement due to fault slip
    ss.ku1(:,k)=u1h(x2GPS,0,y2,ss.y3f(k),Wf(k));
end

% Stress kernels from shear zones
ss.k1312=zeros(ss.Nx*ss.Nz); 
ss.k1313=zeros(ss.Nx*ss.Nz);
ss.k1212=zeros(ss.Nx*ss.Nz);
ss.k1213=zeros(ss.Nx*ss.Nz);

ss.k1212f=zeros(length(ss.y3f),ss.Nx*ss.Nz);
ss.k1312f=zeros(length(ss.y3f),ss.Nx*ss.Nz);
ss.k1213f=zeros(length(ss.y3f),ss.Nx*ss.Nz);
ss.k1313f=zeros(length(ss.y3f),ss.Nx*ss.Nz);

% Displacement kernels
ss.ku12=zeros(length(x2GPS),ss.Nx*ss.Nz);
ss.ku13=zeros(length(x2GPS),ss.Nx*ss.Nz);

% Fields from Shear zones
% Evaluate stress at the center of the fault and shear zone patches
% and displacement at the top of the patch
for kx=1:ss.Nx
    for kz=1:ss.Nz
        % Stress on shear zones due to strain in shear zones
        ss.k1212(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1312(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1213(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1313(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));

        % Stress on faults due to strain in shear zones
        ss.k1212f(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),y2-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1312f(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),y2-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1213f(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),y2-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1313f(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),y2-ss.x2c(kx),ss.y3f(:)+dz/2);
        
        % Displacement due to strain in shear zones
        ss.ku12(:,(kz-1)*ss.Nx+kx)=uk12(ss.polesz(kz),L(kx),W(kz),x2GPS(:)-ss.x2c(kx),0*x2GPS(:));
        ss.ku13(:,(kz-1)*ss.Nx+kx)=uk13(ss.polesz(kz),L(kx),W(kz),x2GPS(:)-ss.x2c(kx),0*x2GPS(:));
    end
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
ss.b = ss.a+2.1e-4*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.012*ones(size(ss.y3f));

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
top    = floor(5e3/(Transition/ss.M));
bottom = ceil(15e3/(Transition/ss.M));
ss.b(1:top)      = ss.a(1:top)     -2.1e-4*ones(top,1);
ss.b(bottom:end) = ss.a(bottom:end)-2.1e-4*ones(length(ss.a(bottom:end)),1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Values used for wet olivine from Hirth, G. and D. Kohlstedt (2003)

% Driving strain rate (1/s)
ss.e12p_plate = 1e-14*ones(ss.Nx*ss.Nz,1);
ss.e13p_plate =      zeros(ss.Nx*ss.Nz,1);

% Rheological Parameters for diffusion and dislocation creep
% Reference Strain Rate (for stress in MPa)
ss.Adif = 1e6*ones(ss.Nz*ss.Nx,1);
ss.Adis = 90 *ones(ss.Nz*ss.Nx,1);

% Power-Law Exponent for dislocation creep
ss.n = 3.5*ones(ss.Nz*ss.Nx,1);

% Activation Energy Wet Oliving (J/mol)
ss.Qdif = 335e3*ones(ss.Nz*ss.Nx,1);
ss.Qdis = 480e3*ones(ss.Nz*ss.Nx,1);

% Activation Volume (m^3/mol)
ss.Voldif = 4e-6*ones(ss.Nz*ss.Nx,1);  
ss.Voldis = 11e-6*ones(ss.Nz*ss.Nx,1);  

% Grain size (m) for diffusion creep
ss.d    = 1e-2*ones(ss.Nz*ss.Nx,1);
ss.pexp = 3*ones(ss.Nz*ss.Nx,1);

% Water fugacity (H/10^6 Si)
ss.COH = 1000*ones(ss.Nz*ss.Nx,1);
ss.r   = 1.2*ones(ss.Nz*ss.Nx,1);

% Temperature (K) and Confining pressure (Pa)
k       = 3.138;         % thermal conductivity (W/m/K)
Cp      = 1171 ;         % specific heat (J/kg/K)
Rm      = 3330 ;         % mantle density (kg/m^3)
Kappa   = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
t_plate = 2e15;          % Plate age (s)

ss.Tprof  = 300+1380*erf(ss.x3c/(sqrt(4* Kappa * t_plate)));  % Kelvin
Te0 = repmat(ss.Tprof',ss.Nx,1);
Te0 = reshape(Te0,[ss.Nx*ss.Nz,1]);

Pconf = Rm*9.8*ss.x3c;  
ss.P  = repmat(Pconf',ss.Nx,1);
ss.P  = reshape(ss.P,[ss.Nx*ss.Nz,1]);

% Coefficients for dislocation and diffusion creep
ss.Const_dis  = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)).*ss.COH.^(ss.r);
ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)).*ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

% Strengh profile
s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
s130 = zeros(size(s120));
e120 = zeros(size(s120));
e130 = zeros(size(s120));

% Fault Strength
ss.strength = ss.sigmab.*(ss.mu0+(ss.a-ss.b).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);

% Plot strength profiles
if false
    figure(1);clf;
    subplot(2,1,1);
    plot(ss.y3f/1e3,ss.strength)
    xlabel('Depth (km)')
    ylabel('Strength (MPa)');
    subplot(2,1,2);
    plot(ss.x3c/1e3,log10(s120(1:ss.Nx:end)))
    xlim([ss.x3c(1)/1e3,ss.x3c(end)/1e3]);
    xlabel('Depth (km)')
    ylabel('Strength (MPa) log10');
    return
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% state parameters
ss.dgfF=4;
ss.dgfS=4;

%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF+ss.Nx*ss.Nz*ss.dgfS,1);

% Fault patches
Y0(1:ss.dgfF:ss.M*ss.dgfF) = zeros(size(ss.y3f));
Y0(2:ss.dgfF:ss.M*ss.dgfF) = ss.strength;
Y0(3:ss.dgfF:ss.M*ss.dgfF) = log(ss.Vo./ss.V_plate);
Y0(4:ss.dgfF:ss.M*ss.dgfF) = log(ss.V_plate*0.98./ss.Vo);

% Shear zones
Y0(ss.M*ss.dgfF+1:ss.dgfS:end) = s120;
Y0(ss.M*ss.dgfF+2:ss.dgfS:end) = s130;
Y0(ss.M*ss.dgfF+3:ss.dgfS:end) = e120;
Y0(ss.M*ss.dgfF+4:ss.dgfS:end) = e130;

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) odeViscoelastic_acceleration(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6); 
[t,Y]=ode45(yp,[0 1e10],Y0,options);
toc
%%
% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end

% Strain rate at center (x2=0)
Ep=sqrt(Yp(:,ss.M*ss.dgfF+floor(ss.Nx/2)*ss.dgfS+3:ss.dgfS*ss.Nx:end)'.^2 +...
        Yp(:,ss.M*ss.dgfF+floor(ss.Nx/2)*ss.dgfS+4:ss.dgfS*ss.Nx:end)'.^2);
    
% Slip Velocity
V=repmat(ss.Vo',size(Y,1)-1,1).*exp(Y(2:end,4:ss.dgfF:ss.M*ss.dgfF));

% Maximum Velocity
Vmax=max(V,[],2);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
f2 = figure(2);clf;set(gcf,'name','Time Evolution')
f2a = subplot(4,1,1);cla;
pcolor(t(1:end-1)/3.15e7,ss.y3f/1e3,log10(V')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
colormap(f2a,parula);
title(h,'Slip Rate (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

f2b = subplot(4,1,2);cla;
pcolor(t(1:end-1)/3.15e7,ss.x3c/1e3,log10(Ep)), shading flat
set(gca,'YDir','reverse');

caxis([log10(min(min(Ep-realmin))) log10(max(max(Ep+realmin)))]);
h1=colorbar('Location','NorthOutside');
colormap(f2b,hot);
title(h1,'Strain Rate (1/s)')
xlabel('Time (Yr)')
ylabel('Depth (km)');

f2c = subplot(4,1,3);cla;
plot(t(1:end-1)/3.15e7,log10(Vmax))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Maximum slip rate on fault')

f2d = subplot(4,1,4);cla;
plot(t(1:end-1)/3.15e7,log10(V(:,floor((top+bottom)/2))))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zone')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
f3 = figure(3);clf;set(gcf,'name','Time Step Evolution')

f3a = subplot(4,1,1);cla;
pcolor(1:length(t)-1,ss.y3f/1e3,log10(V')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V-realmin))) max(max(log10(V+realmin)))]);
colormap(f3a,parula);
title(h,'Slip Rate (m/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

f3b = subplot(4,1,2);cla;
pcolor(1:length(t)-1,ss.x3c/1e3,log10(Ep)), shading flat
set(gca,'YDir','reverse');

caxis([log10(min(min(Ep-realmin))) log10(max(max(Ep+realmin)))]);
h1=colorbar('Location','NorthOutside');
colormap(f3b,hot);
title(h1,'Strain Rate (1/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

f3c = subplot(4,1,3);cla;
plot(1:length(t)-1,log10(Vmax))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Maximum slip rates on fault')

f3d = subplot(4,1,4);cla;
plot(1:length(t)-1,log10(V(:,floor((top+bottom)/2))))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zone')
