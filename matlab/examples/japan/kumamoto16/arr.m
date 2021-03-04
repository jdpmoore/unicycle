function out=arr(T)
R=8.314;
Q=250e3;
out=exp(Q./(R*T))/(365*24*60*60)/1e6;
end %end arrhenius generator