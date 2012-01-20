function z0 = compute_z0(x) 

global a0;
global a1;
global lambda;
z0 = a0 + a1*cos(2*pi*(x-lambda/2)/lambda);
