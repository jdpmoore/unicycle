
% prescribe strain in the shear-zone-centric coordinate system
epsv11=1e-3;
epsv12=0e-3;
epsv13=0e-6;
epsv22=-1e-3;
epsv23=0e-6;
epsv33=0e-3;

% dimension
L=2; % width down dip
W=2; % width down dip
T=2; % thickness

% position
q1=-1;
q2=0;
q3=1; % depth

% orientation
theta=0;

% elastic structure
G=1;     % rigidity
nu=0.25; % Poisson's ratio

% observation points
N=256;
dx1=0.05;
dx2=0.05;
x1=repmat((-N/2:(N/2)-1)'*dx1,1,N);
x2=repmat((-N/2:(N/2)-1)*dx2,N,1);
x3=0*x2+0;

tic
%[u1n2,u2n2,u3n2]=computeDisplacementVerticalShearZoneSurfaceTanhSinh( ...
%    x1,x2,   q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);

%[u1n,u2n,u3n]=computeDisplacementVerticalShearZoneTanhSinh( ...
%    x1,x2,x3,q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);

[u1a,u2a,u3a]=unicycle.greens.computeDisplacementVerticalShearZone( ...
    x1,x2,x3,q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);
toc

u1n=u1a;
u2n=u2a;
u3n=u3a;

%% plot displacements

figure(1);clf;
subplot(3,2,1);cla;hold on;
pcolor(x2/1e3,x1/1e3,u1a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim1=get(gca,'clim');
ylabel(h,'u_1');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution u1');

subplot(3,2,3);cla;hold on;
pcolor(x2/1e3,x1/1e3,u2a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim2=get(gca,'clim');
ylabel(h,'u_2');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution u2');


subplot(3,2,5);cla;hold on;
pcolor(x2/1e3,x1/1e3,-u3a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim3=get(gca,'clim');
ylabel(h,'u_3');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution u3');

subplot(3,2,2);cla;hold on;
pcolor(x2/1e3,x1/1e3,u1n-u1a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'u_1');
set(gca,'clim',clim1/1e14);
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals u1');

subplot(3,2,4);cla;hold on;
pcolor(x2/1e3,x1/1e3,u2n-u2a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'u_2');
set(gca,'clim',clim2/1e14)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals u2');


subplot(3,2,6);cla;hold on;
pcolor(x2/1e3,x1/1e3,-u3n+u3a),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'u_3');
set(gca,'clim',clim3/1e14)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals u3');

%%
figure(62);clf;
cla;hold on;
axis equal tight
down=4;
set(gca,'Clipping','on');
pcolor(x2/1e3,x1/1e3,-u3a),shading flat;
quiver(x2(1:down:end,1:down:end)/1e3,x1(1:down:end,1:down:end)/1e3, ...
    u2a(1:down:end,1:down:end),u1a(1:down:end,1:down:end));
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'Displacement');
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim'))));
%set(gca,'xlim',[-1 1]*60,'ylim',[0 1]*120);
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
