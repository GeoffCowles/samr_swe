function [ierr] = pencilplot(casename,iskip);
close all
%casename = 'roelvink';
dir = ['../pencil_' casename];  
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(4)/2 scrsz(3)/2]) 
C_manning = .01;

%if(exist(option));
% if(option == -1)
%   last = true;
% end;
%end;
eps = 1e-5;
if(~exist('iskip'))
 iskip = 1;
end;


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
  [x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
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
  end;

  clf
  drange = max(max(dinit)-min(dinit),max(dinit));
  subplot(7,1,1)
  plot(xinit,dinit,'r'); hold on
  plot(x,d,'k-+')
  title('d')
  axis([min(x),max(x),min(d)-.1*drange,max(d)+.1*drange]);
  subplot(7,1,2)
  plot(xinit,uinit,'r'); hold on
  plot(x,u,'k')
  axis([min(x),max(x),min(u)-.1,max(u)+.1]);
  title('u')
  subplot(7,1,3)
  plot(xinit,etainit,'r'); hold on
  plot(x,eta,'k')
  axis([min(x),max(x),.9*min(eta-eps),1.1*max(eta+eps)]); 
  title('eta')
  subplot(7,1,4)
  plot(xinit,binit,'r'); hold on
  plot(x,b,'k')
  plot(x,b+d,'r')
  axis([min(x),max(x),min(b)-.1,max(b+d)+.1]); 
  title('b')
  subplot(7,1,5)
  plot(xinit,wdinit,'r'); hold on
  plot(x,wd,'k')
  title('wd')
  axis([min(x),max(x),0,1.2]);
  subplot(7,1,6)
  dx = diff(x);
  plot(xinit,[dxinit',dxinit(end)'],'r'); hold on
  plot(x,[dx',dx(end)'],'k'); hold on
  axis([min(x),max(x),min(min(dxinit),min(dx))-.01*min(dx),max(max(dxinit),max(dx))+.01*min(dx)]); 
  title('dx')
  subplot(7,1,7)
  plot(x,ST_manning(d,u,C_manning),'k')
  title('taub (N)')
%  axis([min(),max(x),0,1.2]);

  fprintf('%d %f %f\n',i,max(abs(u)),max(eta))
  pause
  else
  error('end of data files\n')
  end;

end;

  
