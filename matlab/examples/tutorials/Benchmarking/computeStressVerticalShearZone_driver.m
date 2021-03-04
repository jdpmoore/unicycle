
% prescribe strain in the shear-zone-centric coordinate system
epsv11=1e-6;
epsv12=0e-6;
epsv13=0e-6;
epsv22=1e-6;
epsv23=0e-6;
epsv33=1e-6;

% dimension
L=5e3; % width down dip
W=5e3; % width down dip
T=5e3; % thickness

% position
q1=-2.5e3;
q2=0e3;
q3=0e3; % depth

% orientation
theta=0;

% elastic structure
G=1;     % rigidity
nu=0.25; % Poisson's ratio

% observation points
N=64;
dx1=5e3;
dx2=5e3;
x1=repmat((-N/2:(N/2)-1)'*dx1,1,N);
x2=repmat((-N/2:(N/2)-1)*dx2,N,1);
x3=0*x2+q3+W/2;

tic
[u1a,u2a,u3a]=computeDisplacementVerticalShearZone( ...
    x1,x2,x3,q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);

[s11n,s12n,s13n,s22n,s23n,s33n]=computeStressVerticalShearZoneFiniteDifference( ...
    x1,x2,x3,q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);

[s11a,s12a,s13a,s22a,s23a,s33a]=computeStressVerticalShearZone( ...
    x1,x2,x3,q1,q2,q3,L,T,W,theta,epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu);
toc

s11pa=(cosd(-theta)*s11a-sind(-theta)*s12a)*cosd(theta)-(cosd(theta)*s12a-sind(-theta)*s22a)*sind(-theta);
s12pa=(cosd(-theta)*s11a-sind(-theta)*s12a)*sind(-theta)+(cosd(theta)*s12a-sind(-theta)*s22a)*cosd(theta);
s13pa= cosd(-theta)*s13a-sind(-theta)*s23a;
s22pa=(sind(-theta)*s11a+cosd(theta)*s12a)*sind(-theta)+(sind(-theta)*s12a+cosd(theta)*s22a)*cosd(theta);
s23pa= sind(-theta)*s13a+cosd(theta)*s23a;
s33pa=s33a;

s11pn=(cosd(-theta)*s11n-sind(-theta)*s12n)*cosd(theta)-(cosd(theta)*s12n-sind(-theta)*s22n)*sind(-theta);
s12pn=(cosd(-theta)*s11n-sind(-theta)*s12n)*sind(-theta)+(cosd(theta)*s12n-sind(-theta)*s22n)*cosd(theta);
s13pn= cosd(-theta)*s13n-sind(-theta)*s23n;
s22pn=(sind(-theta)*s11n+cosd(theta)*s12n)*sind(-theta)+(sind(-theta)*s12n+cosd(theta)*s22n)*cosd(theta);
s23pn= sind(-theta)*s13n+cosd(theta)*s23n;
s33pn=s33n;


% s11pn=s11n;
% s12pn=s12n;
% s13pn=s13n;
% s22pn=s22n;
% s23pn=s23n;
% s33pn=s33n;
% 
% s11pa=s11a;
% s12pa=s12a;
% s13pa=s13a;
% s22pa=s22a;
% s23pa=s23a;
% s33pa=s33a;

%% plot displacements

figure(1);clf;
subplot(3,4,1);cla;hold on;
pcolor(x2/1e3,x1/1e3,s11pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim11=get(gca,'clim');
ylabel(h,'s_{11}');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s11');

subplot(3,4,3);cla;hold on;
pcolor(x2/1e3,x1/1e3,s12pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim12=get(gca,'clim');
ylabel(h,'s_{12}');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s12');


subplot(3,4,5);cla;hold on;
pcolor(x2/1e3,x1/1e3,s13pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim13=get(gca,'clim');
ylabel(h,'s_{13}');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s13');

subplot(3,4,7);cla;hold on;
pcolor(x2/1e3,x1/1e3,s22pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim22=get(gca,'clim');
ylabel(h,'s_{22}');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s22');

subplot(3,4,9);cla;hold on;
pcolor(x2/1e3,x1/1e3,s23pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{23}');
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim23=get(gca,'clim');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s23');

subplot(3,4,11);cla;hold on;
pcolor(x2/1e3,x1/1e3,s33pa+s22pa+s11pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{33}');
set(gca,'clim',[-1 1]*max(abs(get(gca,'clim')))/1e0);
clim33=get(gca,'clim');
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Analytic solution s33');


scale=1e6;

subplot(3,4,2);cla;hold on;
pcolor(x2/1e3,x1/1e3,s11pn-s11pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'u_2');
set(gca,'clim',clim11/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s11');

subplot(3,4,4);cla;hold on;
pcolor(x2/1e3,x1/1e3,s12pn-s12pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'u_2');
set(gca,'clim',clim12/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s12');

subplot(3,4,6);cla;hold on;
pcolor(x2/1e3,x1/1e3,s13pn-s13pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{13}');
set(gca,'clim',clim13/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s13');

subplot(3,4,8);cla;hold on;
pcolor(x2/1e3,x1/1e3,s22pn-s22pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{22}');
set(gca,'clim',clim22/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s22');

subplot(3,4,10);cla;hold on;
pcolor(x2/1e3,x1/1e3,s23pn-s23pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{23}');
set(gca,'clim',clim23/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s23');

subplot(3,4,12);cla;hold on;
pcolor(x2/1e3,x1/1e3,s33pn-s33pa),shading flat;
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'s_{33}');
set(gca,'clim',clim33/scale)
axis equal tight
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
title('Residuals s33');


%%
figure(62);clf;
cla;hold on;
axis equal tight
down=2;
set(gca,'Clipping','on');
pcolor(x2/1e3,x1/1e3,-u3a),shading flat;
quiver(x2(1:down:end,1:down:end)/1e3,x1(1:down:end,1:down:end)/1e3, ...
    u2a(1:down:end,1:down:end),u1a(1:down:end,1:down:end));
plot(q2/1e3,q1/1e3,'+')
h=colorbar();
ylabel(h,'Displacement');
%set(gca,'xlim',[-1 1]*60,'ylim',[0 1]*120);
box on, grid on
xlabel('x_2 (km)');
ylabel('x_1 (km)');
