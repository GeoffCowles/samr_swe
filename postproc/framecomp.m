function ierr = framecomp(case1,num1,case2,num2)
close all
dir1 = ['../pencil_' case1]; 
dir2 = ['../pencil_' case2]; 
%label1 = '1st order';
%label2 = 'minmod';
label1 = case1;
label2 = case2 ;

i = num1;
if(i < 10)
  istr = ['000' int2str(i)];
elseif(i < 100) & (i >= 10)
  istr = ['00' int2str(i)];
elseif(i < 1000) & (i >= 100)
  istr = ['0' int2str(i)];
else
  istr = [int2str(i)];
end;
i = num2;
if(i < 10)
  istr2 = ['000' int2str(i)];
elseif(i < 100) & (i >= 10)
  istr2 = ['00' int2str(i)];
elseif(i < 1000) & (i >= 100)
  istr2 = ['0' int2str(i)];
else
  istr2 = [int2str(i)];
end;

 
fname = ['./' dir1 '/' 'pencil_' istr '.dat'];
if(~exist(fname))
  error(['fname: ' fname ' does not exist']);
end;
[x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
if(~exist(fname))
  error(['fname: ' fname ' does not exist']);
end;
fname = ['./' dir2 '/' 'pencil_' istr2 '.dat'];
[x2,d2,u2,eta2,b2,wd2] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);

figure
plot(x,d,'r+'); hold on
plot(x2,d2,'k')
legend(label1,label2)
title('d')

figure
plot(x,u,'r+'); hold on
plot(x2,u2,'k')
legend(label1,label2)
title('u')

