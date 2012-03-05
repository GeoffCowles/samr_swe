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


t = 3.0;
Msize = 3.;
simline = 'k';%'k-o';
iskip = 1;
iskip2 = 8;
Lwidth = 1.2;
z = zeros(numel(xg),1);
fname = 'exner_3.0s.dat';
[xi,di,ui,etai,b,wd,blval] = textread(fname ,'%f %f %f %f %f %f %f\n','headerlines',1);
xi = xi(1:iskip:end);
b  = b (1:iskip:end);
for i=1:numel(xg);
 x = xg(i);
 z(i) = fzero('exner',0) ;
end;
subplot(2,1,1)
plot(xg,compute_z0(xg),'b--','LineWidth',Lwidth); hold on;
plot(xg,z,'r','LineWidth',Lwidth);
plot([-1000,-998],[-1000,-998],'k-o');
legend('Initial Condition','Analytical Solution','Computed Solution')

plot(xi,b,simline,'LineWidth',Lwidth,'MarkerSize',Msize); 
plot(xi(1:iskip2:end),b(1:iskip2:end),'ko','LineWidth',Lwidth,'MarkerSize',Msize); 
axis([0,20,0,3.0])
xlabel('x(m)','FontSize',14)
ylabel('z_b(m)','FontSize',14)


t = 4.5;
z = zeros(numel(xg),1);
fname = 'exner_4.5s.dat';
[xi,di,ui,etai,b,wd,blval] = textread(fname ,'%f %f %f %f %f %f %f\n','headerlines',1);
xi = xi(1:iskip:end);
b  = b (1:iskip:end);
for i=1:numel(xg);
 x = xg(i);
 z(i) = fzero('exner',0) ;
end;
subplot(2,1,2)
plot(xg,compute_z0(xg),'b--','LineWidth',Lwidth); hold on;
plot(xg,z,'r','LineWidth',Lwidth);
plot(xi,b,simline,'LineWidth',Lwidth,'MarkerSize',Msize); 
plot(xi(1:iskip2:end),b(1:iskip2:end),'ko','LineWidth',Lwidth,'MarkerSize',Msize); 
axis([0,20,0,3.0])
xlabel('x(m)','FontSize',14)
ylabel('z_b(m)','FontSize',14)
