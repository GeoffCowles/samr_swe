% compute analytical solution to the exner equation as a function of time
% see kubatko et al, ocean modelling 2006
clear all; close all;

xg = 0:.05:20;
global a0; a0 = 1;
global a1; a1 = 1; 
global lambda; lambda = 20;
global zeta; zeta = 3;
global Aqf; Aqf = 1;
global x
global t


t = 3.0
z = zeros(numel(xg),1);
fname = 'exner_3.0s.dat'
[xi,di,ui,etai,b,wd,blval] = textread(fname ,'%f %f %f %f %f %f %f\n','headerlines',1);
xi = xi(1:3:end);
b  = b (1:3:end);
for i=1:numel(xg);
 x = xg(i);
 z(i) = fzero('exner',0) ;
end;
subplot(2,1,1)
plot(xg,compute_z0(xg),'b--',xg,z,'b',xi,b,'ko'); 
axis([0,20,0,2.5])


t = 4.5
z = zeros(numel(xg),1);
fname = 'exner_4.5s.dat'
[xi,di,ui,etai,b,wd,blval] = textread(fname ,'%f %f %f %f %f %f %f\n','headerlines',1);
xi = xi(1:3:end);
b  = b (1:3:end);
for i=1:numel(xg);
 x = xg(i);
 z(i) = fzero('exner',0) ;
end;
subplot(2,1,2)
plot(xg,compute_z0(xg),'b--',xg,z,'b',xi,b,'ko'); 
axis([0,20,0,2.5])
