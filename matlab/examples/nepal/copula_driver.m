% This is copula driver to optimize the parameter search.
% Author : Sagar Masuti.
clear all;
close all;

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end
eval(['addpath ' home '/Documents/src/unicycle/matlab'])

import unicycle.*

%addpath ../../../
import unicycle.optim.*

config_filepath='./';
% ----------------------------------1------------------------------------------
% Create instance of copula.
cop=copula(config_filepath,'config.dat');
cop.v=0.0001;

% Read gps network station names from network file. For example, gps_network.dat 
% contains columns as below. 
%               n station north  east  up
%               1    GPS1  -200   100   0 
%               2    GPS2  -160   100   0 

gpsstring='gps/';            
gpsfilename=strcat(gpsstring,cop.gps);
[~,stationName,~,~,~]=textread(gpsfilename,'%d %s %f %f %f','commentstyle','shell');

% ----------------------------------2------------------------------------------
% Load the observed data from the files.
d = [];ts=[];
for i = 1:length(stationName)
    filename=strcat(cop.location,strjoin(stationName(i)),'.ned');
    fid=fopen(filename,'r');
    out=textscan(fid,'%f %f %f %f %f %f %f','commentstyle','#');
    time=out{1};
    north=out{2}; 
    east=out{3};
    up=out{4};
    north_err=out{5};
    east_err=out{6};
    up_err=out{7};
    north=north+north_err;
    east=east+east_err; 
    up=up+up_err;
    ts=[ts;time];       
    d=[d;north];
    d=[d;east];
    d=[d;up];
    fclose(fid);
end
tpts=unique(sort(ts));

% ----------------------------------3------------------------------------------
% Call the copula posterior method to get the posterior mean and covariance.
[mp,Sm]=cop.posterior(d);


% ----------------------------------4------------------------------------------
%% Plotting the results. 
Ns=cop.Ns;
M=cop.M;

prior_dist_type=cop.prior_dist_type;      
scale=@(x) 10.^x;
bounds=cop.bounds;
%%
% If number of parameters are more than 2, then choose id for which you want 
% to plot the results.
id = [1,2].';
S = Sm(id,id);
S = (S+S.')/2;
mu = mp(id).';

m1 = linspace(bounds(id(1),1)+1e-4,bounds(id(1),2)-1e-4,1000);
m2 = linspace(bounds(id(2),1)+1e-4,bounds(id(2),2)-1e-4,1000);
nm1 = length(m1);
nm2 = length(m2);
om1=m1;
om2=m2;
[m1,m2] = meshgrid(m1,m2);

urnd=haltonset(M,'Skip',1e3,'Leap',1e2);
urnd=scramble(urnd,'RR2');
urnd=qrandstream(urnd);
msCDF=qrand(urnd,Ns);
ms=zeros(Ns,M);
for i=1:M
    ms(:,i)=bounds(i,1)+(bounds(i,2)-bounds(i,1))*msCDF(:,i);
end
 
fm1 = 1/(bounds(id(1),2)-bounds(id(1),1))*ones(nm1*nm2,1);
fm2 = 1/(bounds(id(2),2)-bounds(id(2),1))*ones(nm1*nm2,1);
gm1 = norminv((m1(:)-bounds(id(1),1))/(bounds(id(1),2)-bounds(id(1),1)));
gm2 = norminv((m2(:)-bounds(id(2),1))/(bounds(id(2),2)-bounds(id(2),1)));

mpdf = mvnpdf([gm1,gm2],mu,S)./normpdf(gm1)./normpdf(gm2).*fm1.*fm2;
mpdfre=reshape(mpdf,nm1,nm2);

[~,posx]=max(max(mpdfre));
[~,posy]=max(max(mpdfre,[],2));

fig=figure(2);clf; 
set(gcf,'Color','w');
mainAx = axes('Parent', fig, 'Position', [0.1, 0.1, 0.5, 0.5]);
pcolor(scale(m1),scale(m2),mpdfre);   % posterior distribution
hold on
scatter(scale(ms(:,id(1))),scale(ms(:,id(2))),10,'white');
plot(scale([om1(1) om1(end)]),scale(om2(posy))*[1 1],'r-.')
plot(scale(om1(posx))*[1 1],scale([om2(1) om2(end)]),'r-.')
shading interp
xlabel('Frictional parameter (a)');
ylabel('Effective Viscosity (Pa s)');

margeAx1 = axes('Parent', fig, 'Position', [0.1,0.63,0.5,0.25]);
marginalm1=trapz(om2,mpdfre,1);
nc=trapz(scale(om1),marginalm1,2);
marginalm1=marginalm1/nc;

m1c=mpdfre(posy,:);
plot (scale(om1),marginalm1,'b', 'LineWidth', 2, 'Parent', margeAx1);
set(gca, 'XTick', []);
axis tight;
ylabel('Probability density');

margeAx2 = axes('Parent', fig, 'Position', [0.63,0.1,0.20,0.5]);
marginalm2=trapz(om1,mpdfre,2);
nc=trapz(scale(om2),marginalm2,1);
marginalm2=marginalm2/nc;

m2c=mpdfre(:,posx);
plot (marginalm2,scale(om2),'b', 'LineWidth', 2, 'Parent', margeAx2);
ylim(mainAx.YLim);
set(gca, 'YTick', []);
xlabel('Probability density');
