function [xc, km, kp] = flop(xc)
x = sh_surface(xc);
[x, km, kp] = flop(x);
km = logical(km);
kp = logical(kp);
xc = x.xc;