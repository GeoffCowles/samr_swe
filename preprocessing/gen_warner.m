% generate the patch coordinates for a Roelvink case
clear all;
n = 8; %cells/km
full = 0;

dx = 1./n;
dy = 1./n; %must be integer fraction 1
wallthick = 1;
inlet = 2;
lx = 15.;
ly = 14.;


if(full == 1)
block1 = [0,0,7.5/dx-1,6./dy-1];
block2 = [0,6./dy-1+1,7.5/dx-1+wallthick,7.5/dy-1+1 + inlet/dy-1];
block3 = [0,6./dy-1+1 + inlet/dy-1 + 1,7.5/dx-1,ly./dy-1];  
block4 = [7.5/dx-1+wallthick+1,0,lx/dx-1,ly/dy-1];

block1 = round(block1);
block2 = round(block2);
block3 = round(block3);
block4 = round(block4);

fprintf('domain_boxes = [ (%d,%d) , (%d,%d) ],\n',block1(1),block1(2),block1(3),block1(4));
fprintf('               [ (%d,%d) , (%d,%d) ],\n',block2(1),block2(2),block2(3),block2(4));
fprintf('               [ (%d,%d) , (%d,%d) ],\n',block3(1),block3(2),block3(3),block3(4));
fprintf('               [ (%d,%d) , (%d,%d) ]\n',block4(1),block4(2),block4(3),block4(4));
else;

block1 = [0,0,7.5/dx-1+wallthick,(inlet/2)/dy-1];  
block2 = [0,(inlet/2)/dy-1+1,7.5/dx-1,(ly/2)/dy-1];
block3 = [7.5/dx-1+wallthick+1,0,(lx)/dx-1,(ly/2)/dy-1];

block1 = round(block1);
block2 = round(block2);
block3 = round(block3);

fprintf('domain_boxes = [ (%d,%d) , (%d,%d) ],\n',block1(1),block1(2),block1(3),block1(4));
fprintf('               [ (%d,%d) , (%d,%d) ],\n',block2(1),block2(2),block2(3),block2(4));
fprintf('               [ (%d,%d) , (%d,%d) ] \n',block3(1),block3(2),block3(3),block3(4));

end;



