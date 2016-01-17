function [t,p,r] = kk_cart2sph(u, v, w)

[theta, phi,r] = cart2sph(u, v, w);
p = theta;
t = pi/2-phi;