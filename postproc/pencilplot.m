function [ierr] = pencilplot(casename);
%casename = 'nomotion';
close all
dir = ['../pencil_' casename];  

%if(exist(option));
% if(option == -1)
%   last = true;
% end;
%end;


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
  [x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
  if(i==0)
    dinit = d;
    uinit = u;
    etainit = eta;
    binit =b;
    wdinit = wd;
    dxinit = diff(x);
  end;

  clf
  drange = max(max(dinit)-min(dinit),max(dinit));
  subplot(6,1,1)
  plot(x,dinit,'r'); hold on
  plot(x,d,'k-+')
  title('d')
  axis([min(x),max(x),min(d)-.1*drange,max(d)+.1*drange]);
  subplot(6,1,2)
  plot(x,uinit,'r'); hold on
  plot(x,u,'k')
  axis([min(x),max(x),min(u)-.1,max(u)+.1]);
  title('u')
  subplot(6,1,3)
  plot(x,etainit,'r'); hold on
  plot(x,eta,'k')
  axis([min(x),max(x),.9*min(eta),1.1*max(eta)]); 
  title('eta')
  subplot(6,1,4)
  plot(x,binit,'r'); hold on
  plot(x,b,'k')
  title('b')
  subplot(6,1,5)
  plot(x,wdinit,'r'); hold on
  plot(x,wd,'k')
  title('wd')
  axis([min(x),max(x),0,1.2]);
  subplot(6,1,6)
  dx = diff(x);
  plot(x,[dxinit',dxinit(end)'],'r'); hold on
  plot(x,[dx',dx(end)'],'k'); hold on
  title('dx')
  fprintf('%d %f %f\n',i,max(abs(u)),max(eta))
  pause
  end;

end;

  
