function [xc, km, kp] = flip(xc)
x = sh_surface(xc);
[x, km, kp] = flip(x);
km = logical(km);
kp = logical(kp);
xc = x.xc;