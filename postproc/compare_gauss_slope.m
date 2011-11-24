close all;
fname = {'gauss_centerslope.dat','gauss_edgeslope.dat','gauss_noslope.dat','gauss_centerslope_sedangle1deg.dat'}; 
plotty = {'r-+','k-o','b-d','g-d'}
for i=1:4
[x,d,u,eta,b,wd,blval] = textread(char(fname(i)),'%f %f %f %f %f %f %f\n','headerlines',1);
plot(x,b,char(plotty(i))); hold on;
end;
legend('center slope','edge slope','no slope','center slope with sed angle at 1 deg'); 
