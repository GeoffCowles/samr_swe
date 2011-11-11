% sampson analytical solution for eta
function [bath,eta] = sampson(x,t);
a = 3000.;
h0 = 10.;
tau = .001;
B = 5.;
g = 9.8016;
p = sqrt(8*g*h0)/a;
s = sqrt(p^2-tau^2)/2.;

% compute bathy at x
bath = h0*((x/a)^2);

% compute eta at (x,t)
fac =  (a^2*B^2*exp(-t*tau))/(8*g^2*h0);
fac2 = (.25*tau^2 - s^2);
fac3 = (x*exp(-t*tau/2))/g;
eta = h0 + fac*(-s*tau*sin(2*s*t) + fac2*cos(2*s*t)) - (B^2*exp(-t*tau))/(4*g) - fac3*(B*s*cos(s*t) + .5*tau*B*sin(s*t));
eta = max(eta,bath);


