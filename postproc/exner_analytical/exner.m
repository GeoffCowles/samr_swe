function resid = exner(z)
  global x;
  global Aqf;
  global zeta;
  global t;
  resid = z-compute_z0(x-Aqf/((zeta-z)^2)*t);

