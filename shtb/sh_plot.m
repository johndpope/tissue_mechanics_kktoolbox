function [xr, yr, zr] = sh_plot(c)
%%% just plot the radial function given its SH coefficients
%%% i.e. we are performing spherical harmonic synthesis
if nargin ==1, gdim = 40;L_max = sqrt(length(c))-1;end
[t wt]                  = gaussquad(gdim, 0, pi);
[p wp]                  = gaussquad(gdim,0,2*pi);
[p t]                   = meshgrid(p,t);
[wp wt]                 = meshgrid(wp, wt);
[Y_LK P_LK]				= precalc_ylk_cos_sin_old(p, t, L_max);
c = c(:)';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
r = sum(c.*Y_LK,3);
xr = r.*sin(t).*cos(p);
yr = r.*sin(t).*sin(p);
zr = r.*cos(t);
surf(xr,yr,zr);daspect([1 1 1]);

