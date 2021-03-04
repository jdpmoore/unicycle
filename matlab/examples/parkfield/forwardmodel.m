function res = forwardmodel(m,evl,time,options,gps)
% funtion forwardmodel computes residuals between unicycle simulations and
% GPS data based on two model parameters (a-b) and Vo in that order.

% rate-and-state coefficients
evl.rcv.b=evl.rcv.a-m(1);
% reference velocity (m/yr)
evl.rcv.Vo=zeros(size(evl.rcv.a))+m(2);

% initial condition
y0=zeros(1,evl.rcv.N*evl.dgf);
opts=unicycle.ode.odeset('Refine',1,'RelTol',1e-15,'InitialStep',1e-5,'Verbose',false);
[evl.t,evl.y]=evl.ode45(time,y0,opts);

res=sqrt(sum((gps.d-gps.dataVectorFromModel(evl)).^2));

end