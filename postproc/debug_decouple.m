function [ierr] = pencilplot(casename,iskip);
close all
%casename = 'roelvink';
dir = ['../pencil_' casename];  
scrsz = get(0,'ScreenSize');
%figure('Position',[1 1 scrsz(4)/2 .66*scrsz(3)]) 
C_manning = .01;
a = 0;

%if(exist(option));
% if(option == -1)
%   last = true;
% end;
%end;
eps = 1e-5;
if(~exist('iskip'))
 iskip = 1;
end;

nplots = 7;


for i=0:iskip:10000
  if(i < 10)
    istr = ['000' int2str(i)]; 
  elseif(i < 100) & (i >= 10)
    istr = ['00' int2str(i)]; 
  elseif(i < 1000) & (i >= 100)
    istr = ['0' int2str(i)]; 
  else(i < 10000) & (i >= 1000)
    istr = [int2str(i)]; 
  end;
 
  fname = ['./' dir '/' 'pencil_' istr '.dat'];
  if(exist(fname))
  [x,d,u,eta,b,wd,blval] = textread(fname,'%f %f %f %f %f %f %f\n','headerlines',1);
  pts = max(numel(find(isnan(d))==1),numel(find(isnan(u))==1)); 
  if(pts > 0)
    error('got nan')
  end;
  if(i==0)
    xinit = x;
    dinit = d;
    uinit = u;
    etainit = eta;
    binit =b;
    wdinit = wd;
    dxinit = diff(x);
    blvalinit = blval;
  end;
  if(max(blval)>.001)&(a==0)
    ; blvalinit = blval; a=1 ; end;
   

  clf
  drange = max(max(dinit)-min(dinit),max(dinit));
%  dx = diff(x);
%  plot(xinit,[dxinit',dxinit(end)'],'r'); hold on
%  plot(x,[dx',dx(end)'],'k'); hold on
%  axis([min(x),max(x),min(min(dxinit),min(dx))-.01*min(dx),max(max(dxinit),max(dx))+.01*min(dx)]);
%  title('dx')

  taux = ST_manning(d,u,C_manning);
  taux = taux-mean(taux); 
  plot(x,taux/max(taux),'r-o'); hold on;
  plot(x,blval/max(blvalinit),'k-+') 
  d = d-mean(d);
  plot(x,d/max(d),'b-d')
  u = u-mean(u);
  plot(x,u/max(u),'g')
  b = b-mean(b);
  plot(x,b/max(binit),'m')
  title('bedlevel (black) and stress (red) and d (blue) and u (green) and b (magenta)')
  axis([.95e4,1.05e4,-1.2,1.2])


  if(nplots >=8)
  subplot(nplots,1,8)
  plot(xinit,wdinit,'r'); hold on
  plot(x,wd,'k')
  title('wd')
  axis([min(x),max(x),0,1.2]);
  end;

%  axis([min(),max(x),0,1.2]);

  fprintf('%d %f %f\n',i,max(abs(u)),max(eta))
  pause
  else
  error('end of data files\n')
  end;

end;

  
