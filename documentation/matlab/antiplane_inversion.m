
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% We use the following coordinate system for the 
% description of the geometry.
%
%                 N x1 (no deformation in this direction)
%                /
%               / 
%   y1,y2,y3 ->@--------------  x2 (Matlab x)
%              |
%              :
%              | 
%              :
%              Z x3 (Matlab y)
%
% Fault slip occurs along the x3 direction. Positive
% displacement is pointing in the x1 direction.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            P H Y S I C A L   M O D E L               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

u1h=@(x2,x3,y2,y3,W) ...
    (-atan2((x2-y2),(x3-y3))+atan2((x2-y2),(x3+y3))...
     +atan2((x2-y2),(x3-y3-W))-atan2((x2-y2),(x3+y3+W)) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            S Y N T H E T I C   D A T A               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% locking depth
L=1;

N=128;
dx2=0.25;
x2=(-N/2:N/2-1)'*dx2;

rng('default');
sigmad=0.01;
d=u1h(x2,0,0,0,L)+2*(rand(N,1)-0.5)*sigmad;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%          G R E E N ' S   F U N C T I O N S           %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

u1h=@(x2,x3,y2,y3,W) ...
    (-atan2((x2-y2),(x3-y3))+atan2((x2-y2),(x3+y3))...
     +atan2((x2-y2),(x3-y3-W))-atan2((x2-y2),(x3+y3+W)) ...
    )/2/pi;

M=20;
dz=0.2;
y3=(0:M-1)'*dz;
W=ones(M,1)*dz;

% forward model
G=zeros(N,M);
for k=1:M
    G(:,k)=u1h(x2,0,0,y3(k),W(k));
end

m1=ones(M,1);
m1(y3>=L)=0;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%       S M O O T H I N G   C O N S T R A I N T S      %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% smoothing constraints
L=diag(ones(M-1,1),1)-diag(ones(M,1),0);
L(M,:)=0;

% smoothing weight
lambda=1e-3; 

% best fit
H=[G;lambda*L]; % combine Green functions
h=[d;zeros(M,1)]; % combine data and smoothing constraint


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                  I N V E R S I O N                   %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% least-squares inversion with positivity constraint
m2=lsqnonneg(G,d);

% least-squares inversion with smoothing and positivity constraint
m3=lsqnonneg(H,h);

% least-squares inversion, no smoothing
m4=(G.'*G)\(G.'*d);

% least-squares inversion, smoothing
m5=(H.'*H)\(H.'*h);

% pseudo-inverse 
tolerance=1e-1;
m6=pinv(G,tolerance)*d;

% least-squares inversion, no smoothing, fewest possible non-zero
% components
m7=H\h;

% pseudo-inverse with smoothing
%tolerance=1e-1;
m8=pinv(H,tolerance)*h;

% sparse Basis pursuit denoising solution
sigma=0.09; % ||Ax - b||_2 < sigma
m9=spgl1.bpdn(G,d,sigma,opts);

% sparse Lasso solution
tau=5; % || m ||_1 < tau
m10=spgl1.lasso(G,d,tau,opts);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                        %
% F O R W A R D   M O D E L   A N D   R E S I D U A L S  %
%                                                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% model predictions
dt=G*m10;

% residuals between observations and best-fitting model
r=d-dt;

% variance reduction
theta=1-sum(r.^2)/sum(d.^2);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                        %
%                     F I G U R E S                      %
%                                                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% plot forward models (the data is a linear combinations of these modes)
figure(1);clf;
subplot(2,1,1);cla
hold on, box on
plot(x2,d,'k-');
plot(x2,G)
set(gca,'ylim',[-1 1]*0.5)
title('Basis functions')
ylabel('Surface displacements')
xlabel('Across-fault distance (km)')
axis tight
legend('Data','basis functions')

subplot(2,1,2);cla
hold on, box on
Lambda=svds(G,M);
plot(log10(Lambda));
plot([1 M],[1 1]*log10(sigmad/1))
axis tight
xlabel('Eigenvalue number')
ylabel('Energy (log10)')
legend('Eigenvalue','SNR')


%% plot data, forward model, residuals

figure(2);clf
subplot(2,1,1);cla
hold on
plot(x2,d,'g+-');
plot(x2,dt,'k-.');
ylabel('Surface displacements')
legend('Data','Best forward model')
title('Data & model')
box on, axis tight

subplot(2,1,2);cla
plot(x2,r,'r');
title(sprintf('Residuals (variance reduction: %6.2f %%)',theta*100))
xlabel('x')
ylabel('Surface displacements')
box on, axis tight

%% target and best-fitting slip distributions
figure(3);clf
hold on
plot(m1,y3,'r')
plot(m3,y3,'k+-')
plot(m10,y3,'g+-.')
legend('Target model','m3','m10','Location','SouthEast')
xlabel('Slip (m)')
ylabel('Depth (km)')
set(gca,'YDir','Reverse','xlim',[-0.2 1.4]);
title('Slip distributions')
box on

