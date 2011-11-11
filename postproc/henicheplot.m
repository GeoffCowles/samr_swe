%function [ierr] = pencilplot(casename);
casename='heniche';
lthck=1.5;
%casename = 'nomotion';
close all
dir = ['../pencil_' casename];  

%if(exist(option));
% if(option == -1)
%   last = true;
% end;
%end;

% load fvcom time
nc = netcdf('fvcom_data/heniche_e.nc');
ftime = nc{'time'}(:);
ftime = ftime*24*3600; %convert to s
x = nc{'x'}(:);
y = nc{'y'}(:);
midpts = find(abs(y-0)<.01);
xpts = x(midpts);


for i=0:10000
  if(i < 10)
    istr = ['000' int2str(i)]; 
  elseif(i < 100) & (i >= 10)
    istr = ['00' int2str(i)]; 
  else(i < 1000) & (i >= 100)
    istr = ['0' int2str(i)]; 
  end;
 
  fname = ['./' dir '/' 'pencil_' istr '.dat'];
  if(exist(fname))
  fid = fopen(fname,'r');
  C = textscan(fid, '%f', 1);
  time = C{1};
  fclose(fid);
  [x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
  tsam(i+1) = time;
  usam(i+1) = u(1);
  if(i==0)
    dinit = d;
    uinit = u;
    etainit = eta;
    binit =b;
  end;
  %find nearest fvcom time
  [diff,frame] = min(abs(ftime-time));
  eta_fvcom = nc{'zeta'}(frame,midpts); 

  clf
  plot(x,b,'r','LineWidth',lthck); hold on
  plot(x,eta,'k','LineWidth',lthck)
  plot(xpts,eta_fvcom+1,'b','LineWidth',lthck);
  axis([min(x),max(x),0,2.0]); 
  title(['T = ' num2str(ceil(time/60.)) ' m'],'FontSize',14)
  %h_legend = legend('bathy','swe-hlle','FVCOM')
  %set(h_legend,'FontSize',14);
  xlabel('x (m)','FontSize',16)
  ylabel('zeta (m)','FontSize',16)

  fprintf('time in minutes %f\n',time/60.);
  pause
  end;

end;

  
figure
ua = nc{'ua'}(:,1);
plot(ftime,ua,'r',tsam,usam,'k');
