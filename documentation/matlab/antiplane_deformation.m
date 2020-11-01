%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Antiplane deformation: displacements and stress.
%
% We use the following coordinate system for the 
% description of the geometry.
%
%                 N x1 (displacement in this direction only)
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

G=1;

s12h=@(x2,x3,y2,y3,W) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-W)./((x2-y2).^2+(x3-y3-W).^2)-(x3+y3+W)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,W) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
   -(x2-y2)./((x2-y2).^2+(x3-y3-W).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;


% fault geometry
s=1;
y2=0;
y3=1e3;
y2=0;
W=1e3;

%% surface displacement profile

N=1024;
dx2=0.05e3;
x2=(-N/2:N/2-1)'*dx2;
u1=s*u1h(x2,0,y2,y3,W);

figure(1);
clf;set(gcf,'Name','Antiplane deformation');
subplot(2,1,1);cla;
hold on, box on
plot(x2,u1), shading flat;
title('Surface displacements');
xlabel('Across-fault distance (km)');
ylabel('Displacement (m)');
axis tight

subplot(2,1,2);cla;
hold on, box on
plot(x2([1,end]),[0 0]);
plot([0 0],[y3 y3+W]);
title('Fault geometry');
xlabel('Across-fault distance (km)');
ylabel('Depth (km)');
set(gca,'YDir','Reverse');
set(gca,'ylim',[-2 y3+W+1]);
set(gca,'xlim',x2([1,end]))
legend('Surface','Asperity')


%% cross section of stress

N=1024;
dx3=0.005e3;
x3=(0:N-1)'*dx3;
s12=s*s12h(0,x3,y2,y3,W);

figure(2);
clf;set(gcf,'Name','Stress change & singularities');
subplot(2,1,1);cla;
hold on, box on
plot(x3/1e3,s12,'r'), shading flat;
title('Stress change');
xlabel('Depth (km)');
ylabel('Stress s_{12} (Pa)');
axis tight

subplot(2,1,2);cla;
hold on, box on
plot(x3/1e3,2*s*u1h(0.001,x3,y2,y3,W));
title('Slip distribution');
xlabel('Depth (km)');
ylabel('Slip (m)');
axis tight

%% displacement cross section

M=128;
N=128;
dx2=0.025e3;
dx3=0.025e3;
[x2,x3]=meshgrid((-N/2:N/2-1)'*dx2,(0:M-1)*dx3);
u1=s*u1h(x2,x3,y2,y3,W);
s12=s*s12h(x2,x3,y2,y3,W);
s13=s*s13h(x2,x3,y2,y3,W);

figure(3);
clf;set(gcf,'Name','Antiplane displacements');
hold on, box on
pcolor(x2/1e3,x3/1e3,u1), shading flat;
contour(x2/1e3,x3/1e3,u1,25,'k-.')
h=colorbar();
ylabel(h,'Slip (m)');
xlabel('Across-fault distance (km)');
ylabel('Depth (km)');
set(gca,'YDir','Reverse');
axis equal tight

figure(4);
clf;set(gcf,'Name','Antiplane stress');
subplot(2,1,1);cla;
hold on, box on
pcolor(x2/1e3,x3/1e3,s12), shading flat;
h=colorbar();
ylabel(h,'s_{12} (m)');
xlabel('Across-fault distance (km)');
ylabel('Depth (km)');
set(gca,'YDir','Reverse');
axis equal tight

subplot(2,1,2);cla;
hold on, box on
pcolor(x2/1e3,x3/1e3,s13), shading flat;
h=colorbar();
ylabel(h,'s_{13} (m)');
xlabel('Across-fault distance (km)');
ylabel('Depth (km)');
set(gca,'YDir','Reverse');
axis equal tight

%% cross section of stress
figure(4);



