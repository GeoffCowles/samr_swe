close all; clear all;
dir1 = '../pencil_dambreaky'; 
dir2 = '../../swe_old/pencil_dambreakx'; 

for i=1:1000

if(i < 10)
  istr = ['000' int2str(i)]; 
elseif(i < 100) & (i >= 10)
  istr = ['00' int2str(i)]; 
else(i < 1000) & (i >= 100)
  istr = ['0' int2str(i)]; 
end;
 
fname = ['./' dir1 '/' 'pencil_' istr '.dat'];
[x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);

fname = ['./' dir2 '/' 'pencil_' istr '.dat'];
[x2,d2,u2,eta2,b2] = textread(fname,'%f %f %f %f  %f\n','headerlines',1);

subplot(2,1,1)
plot(x,d,'r-+'); hold on
plot(x2,d2,'k')
subplot(2,1,2)
plot(x,u,'r-+'); hold on
plot(x2,u2,'k')
title(int2str(i));
pause
clf

  
end;
