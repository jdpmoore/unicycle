
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
% Evaluates the slip history on a spring-slider system %
%     under the rate- and state-dependent friction     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% static friction coefficient
ss.mu0=0.1;
% stiffness
ss.k=0.6;
% frictional parameters
ss.a=1e-2;
ss.b=ss.a+1.2e-3;
% normal stress
ss.sigma=50.0;
% characteristic-weakening distance
ss.L=0.1;
% plate velocity
ss.Vpl=1e-9;
% reference slip rate
ss.Vo=1e-6;

fprintf('k=%f, (a-b)*sigma / L=%f\n', ...
    ss.k,(ss.b-ss.a)*ss.sigma/ss.L);
fprintf('recurrence time (a-b)*sigma / k / Vpl=%f yr\n', ...
    (ss.b-ss.a)*ss.sigma/ss.k/ss.Vpl/3.14e7);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         N U M E R I C A L   S O L U T I O N          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% initialize the function handle with 
% set constitutive parameters
yp=@(t,y) odefun(t,y,ss);

% solve the system
options=odeset('Refine',1,'RelTol',2e-14,'InitialStep',1e-10);
[t,Y]=ode23(yp,[0 3.1e10],[ss.Vpl;ss.mu0*ss.sigma;0],options);

Yp=zeros(size(Y));
for i=1:length(t)
    Yp(i,:)=yp(t(i),Y(i,:));
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

figure(1);clf;set(gcf,'name','Dynamic variables evolution')
subplot(3,1,1);cla;hold on
plot(t/3.14e7,Y(:,1))
plot(t/3.14e7,t*ss.Vpl,'r')
ylabel('slip (m)')
box on
subplot(3,1,2);cla
plot(t/3.14e7,Y(:,2))
ylabel('stress (MPa)')
box on
subplot(3,1,3);cla
plot(t/3.14e7,Y(:,3))
ylabel('ln (theta Vo / L)')
box on


%% phase diagram

figure(2);clf;set(gcf,'name','Phase diagrams')
subplot(2,1,1);cla;hold on
plot(log10(Yp(:,1)),Y(:,2))
plot(log10(ss.Vpl)*[1 1],get(gca,'ylim'))
plot(log10(Yp(:,1)),(ss.mu0+(ss.a-ss.b)*log(Yp(:,1)/ss.Vo))*ss.sigma)
xlabel('Velocity (m/s)')
ylabel('Force (N) log10')
box on
axis tight

subplot(2,1,2);cla;hold on
plot(log10(Yp(:,1)),Y(:,3))
xlabel('Velocity (m/s)')
ylabel('State variable (s) log10')
box on
axis tight