% toro exact riemann case test #1, see toro book page 120ish
clear all; close all;
% plot exact and samrai solutions for toro wet test
dir = '../pencil_dambreakx/';
fname_init = [dir 'pencil_0000.dat']; 
fname_final= [dir 'pencil_0892.dat']; 
[xi,di,ui,etai,b,wd,blval] = textread(fname_init ,'%f %f %f %f %f %f %f\n','headerlines',1);
[xf,df,uf,etaf,b,wd,blval] = textread(fname_final,'%f %f %f %f %f %f %f\n','headerlines',1);


fname = fname_final;
if(exist(fname))
  fid = fopen(fname,'r');
  C = textscan(fid, '%f', 1);
  time = C{1};
  fclose(fid);
end;
[xf,df,uf,etaf,b,wd,blval] = textread(fname_final,'%f %f %f %f %f %f %f\n','headerlines',1);
fprintf('samrai time %f\n',time);

fname = './ana_data/toro_wettest.zeta';
[xzt,zt] = textread(fname,'%f %f\n','headerlines',0);

fname = './ana_data/toro_wettest.u';
[xut,ut] = textread(fname,'%f %f\n','headerlines',0);

subplot(2,1,1)
lwdth = 1.2;
plot(xi,di,'b--','LineWidth',lwdth); hold on;
plot(xzt,zt,'r-.','LineWidth',lwdth); hold on;
plot(xf,df,'k','LineWidth',lwdth)
%h_legend = legend('initial','exact','samr-swe');
%set(h_legend,'FontSize',14);

axis([0,50,-.1,1.1])
ylabel('d(m)','FontSize',14)
xlabel('x(m)','FontSize',14)
%exportfig(gcf,'dambreak_d_wet.eps','format','eps','color','rgb')

%figure
subplot(2,1,2)
plot(xi,ui,'b--','LineWidth',lwdth); hold on;
plot(xut,ut,'r-.','LineWidth',lwdth); hold on;
plot(xf,uf,'k','LineWidth',lwdth)
axis([0,50,-.5,5]);
ylabel('u(m/s)','FontSize',14)
xlabel('x(m)','FontSize',14)
h_legend = legend('initial','exact','simulated');
set(h_legend,'FontSize',14);
%exportfig(gcf,'dambreak_u_wet.eps','format','eps','color','rgb')



