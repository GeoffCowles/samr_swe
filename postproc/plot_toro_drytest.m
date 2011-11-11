clear all; close all;
% plot exact and samrai solutions for toro dry test
fname= './samrai_data/toro_wetdry_init.dat';
[xi,di,ui,etai,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);

fname = './samrai_data/toro_wetdry_final.dat';
if(exist(fname))
  fid = fopen(fname,'r');
  C = textscan(fid, '%f', 1);
  time = C{1};
  fclose(fid);
end;
[xf,df,uf,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
fprintf('samrai time %f\n',time);

fname = './ana_data/toro_drytest.zeta';
[xzt,zt] = textread(fname,'%f %f\n','headerlines',0);

fname = './ana_data/toro_drytest.u';
[xut,ut] = textread(fname,'%f %f\n','headerlines',0);

lwdth = 1.2;
plot(xi,di,'b--','LineWidth',lwdth); hold on;
plot(xzt,zt,'r','LineWidth',lwdth); hold on;
plot(xf,df,'k','LineWidth',lwdth)

%legend('init','exact','hlle')
axis([0,50,-.5,1.5])
ylabel('d','FontSize',14)
xlabel('x','FontSize',14)
exportfig(gcf,'dambreak_d_dry.eps','format','eps','color','rgb')

figure
ui = zeros(prod(size(xi),1));
plot(xi,ui,'k--'); hold on;
plot(xut,ut,'r'); hold on;
plot(xf,uf,'k')
%legend('init','exact','hlle')
axis([0,50,-.5,7])
ylabel('u')
xlabel('x')
exportfig(gcf,'dambreak_u_dry.eps','format','eps','color','rgb')



