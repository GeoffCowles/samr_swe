%function [ierr] = pencilplot(casename);
lthck=1.5;
%casename = 'nomotion';
close all
%dir = ['../pencil_' casename];  

% load fvcom time
nc = netcdf('fvcom_data/heniche_e.nc');
ftime = nc{'time'}(:);
ftime = ftime*24*3600; %convert to s
x = nc{'x'}(:);
y = nc{'y'}(:);
midpts = find(abs(y-0)<.01);
xpts = x(midpts);

frame = 1;
eta_fvcom = nc{'zeta'}(frame,midpts); 
h_fvcom = nc{'h'}(midpts); 
h_fvcom = 1-h_fvcom;

plot(xpts,h_fvcom,'r','LineWidth',lthck); hold on
plot(xpts,eta_fvcom+1,'k','LineWidth',lthck)
plot(xpts,eta_fvcom+1,'b','LineWidth',lthck);
plot(440,1.34,'kd','MarkerSize',14,'MarkerFaceColor','k'); %,'SymbolSize',15)
axis([min(x),max(x),0,2.0]); 
title(['T = ' num2str(ceil(0)) ' m'],'FontSize',14)
  h_legend = legend('bathy','swe-amr','FVCOM')
  set(h_legend,'FontSize',14);
  xlabel('x (m)','FontSize',14)
  ylabel('zeta (m)','FontSize',14)
%exportfig(gcf,'examp.eps','format','eps','color','rgb')
