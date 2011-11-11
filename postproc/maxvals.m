function [ierr] = maxvals(casename);
%casename = 'nomotion';
close all
dir = ['../pencil_' casename];  

%if(exist(option));
% if(option == -1)
%   last = true;
% end;
%end;
nofile = 0;
i = 0
while nofile ==0;
  if(i < 10)
    istr = ['000' int2str(i)]; 
  elseif(i < 100) & (i >= 10)
    istr = ['00' int2str(i)]; 
  else(i < 1000) & (i >= 100)
    istr = ['0' int2str(i)]; 
  end;
 
  fname = ['./' dir '/' 'pencil_' istr '.dat'];
  if(exist(fname))
    i = i + 1;
    fid = fopen(fname,'r');
    C = textscan(fid, '%f', 1);
    time(i) = C{1};
    fclose(fid);
    [x,d,u,eta,b,wd] = textread(fname,'%f %f %f %f %f %f\n','headerlines',1);
    umax(i) = max(abs(u));
  else
    nofile = 1;
  end;
end;

plot(time,umax);
 


  
