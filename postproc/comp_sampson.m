% compare analytical-computed values of sampson case
lthick=1.2;
casename='sampson';
close all
figure
dir = ['../pencil_' casename];

nc = netcdf('fvcom_data/samp_0001.nc');
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
  end;

  %find nearest fvcom time
  [diff,frame] = min(abs(ftime-time));
  eta_fvcom = nc{'zeta'}(frame,midpts); 

   
  plot(x,b,'k--','LineWidth',lthick);  hold on;
  plot(x,eta,'b','LineWidth',lthick);
  nx = prod(size(x));
  for i=1:nx
    [bana(i),eta_ana(i)] = sampson(x(i),time);
  end;
  fprintf('time: %f\n',time);
  plot(xpts,eta_fvcom,'k','LineWidth',lthick);
  plot(x,eta_ana,'r','LineWidth',lthick);
  ylabel('T=1500s','FontSize',14)
  %h_legend=legend('bathymetry','swe-amr','FVCOM-2D','analytical');
  %set(h_legend,'FontSize',14);
  pause
  clf

end;
