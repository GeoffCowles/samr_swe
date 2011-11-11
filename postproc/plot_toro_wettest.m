clear all; close all;
% plot exact and samrai solutions for toro wet test
fname= './samrai_data/toro_wet_init.dat';
[xi,di,ui,etai,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);

fname = './samrai_data/toro_wet_final.dat';
if(exist(fname))
  fid = fopen(fname,'r');
  C = textscan(fid, '%f', 1);
  time = C{1};
  fclose(fid);
end;
[xf,df,uf,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
fprintf('samrai time %f\n',time);

fname = './ana_data/toro_wettest.zeta';
[xzt,zt] = textread(fname,'%f %f\n','headerlines',0);

fname = './ana_data/toro_wettest.u';
[xut,ut] = textread(fname,'%f %f\n','headerlines',0);

lwdth = 1.2;
plot(xi,di,'b--','LineWidth',lwdth); hold on;
plot(xzt,zt,'r','LineWidth',lwdth); hold on;
plot(xf,df,'k','LineWidth',lwdth)
h_legend = legend('initial','exact','swe-amr')
set(h_legend,'FontSize',14);

axis([0,50,-.1,1.1])
ylabel('d','FontSize',14)
xlabel('x','FontSize',14)
exportfig(gcf,'dambreak_d_wet.eps','format','eps','color','rgb')

figure
plot(xi,ui,'k--','LineWidth',lwdth); hold on;
plot(xut,ut,'r','LineWidth',lwdth); hold on;
plot(xf,uf,'k','LineWidth',lwdth)
axis([0,50,-.5,4]);
ylabel('u')
xlabel('x')
exportfig(gcf,'dambreak_u_wet.eps','format','eps','color','rgb')



