function [u, v, w] = kk_sph2cart(t,p,r)
%% sph2cart conversion using the kk convension (also used in the calculation of the SH functions etc...)
phi = pi/2-t;
theta = p;
[u, v, w] = sph2cart(theta,phi,r);