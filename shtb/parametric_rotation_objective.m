function res = parametric_rotation_objective(ang,xc,yc,zc,D, nY1, nY2, plot_flag)
%% calculate the rotation equation results for fsolve
%%% D contains two elements. the y axis value and the  x axis value 
%%% attainable for the smallest and second smallest half axes of the first order ellipsoid, when the other two
%%% coordinates are zero.
if nargin<6, plot_flag = 0;end
ang = real(ang);
ang(1) = mod(ang(1), 2*pi);
ang(2) = mod(ang(2), 2*pi);
ang(3) = mod(ang(3), 2*pi);
%% plot
if plot_flag
    disp(ang);
    [xc1] = sh_rot(xc, ang(1), ang(2), ang(3));
    [yc1] = sh_rot(yc, ang(1), ang(2), ang(3));
    [zc1] = sh_rot(zc, ang(1), ang(2), ang(3));
    X_o1 = cs2nocs([xc1(:)' yc1(:)' zc1(:)']);
    dfig(1);clf;plot_shps(X_o1,3,'red');view(3);lighting gouraud;camlight;
end
%% rotate the spherical harmonics according to ang
[xc] = sh_rot(xc(1:4), ang(1), ang(2), ang(3));
[yc] = sh_rot(yc(1:4), ang(1), ang(2), ang(3));
[zc] = sh_rot(zc(1:4), ang(1), ang(2), ang(3));

%% % evaluate the length of the surface vector at theta = 0 (at the minimum
%%% this is aligned with the longest axis of the object's first order
%%% ellipsoid
% t = 0;p = 0;
% X1 = xc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + xc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + xc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Y1 = yc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + yc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + yc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Z1 = zc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + zc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + zc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
X1 = nY1*xc(2:4);
Y1 = nY1*yc(2:4);
Z1 = nY1*zc(2:4);
St0 = sqrt(X1^2 + Y1^2 + Z1^2);   % surface vector length at theta = 0 (i.e. at north pole)
%% evaluate the length of the surface vector at phi = 0 and theta = pi/2
% t = pi/2;p = 0;
% X2 = xc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + xc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + xc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Y2 = yc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + yc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + yc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Z2 = zc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + zc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + zc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);

X2 = nY2*xc(2:4);
Y2 = nY2*yc(2:4);
Z2 = nY2*zc(2:4);

Sp0 = sqrt(X2^2 + Y2^2 + Z2^2);   % surface vector length at phi = 0

%%% a minimum is achieved if the Northpole (i.e. the point at theta = 0
%%% (and phi = 0)) is along the longest axis of the first-order ellipsoid
%%% and simultaneously the point where the Greenwich meridian meets the
%%% equator (i.e. phi = 0 and theta = pi/2), is along the shortest axis.

res = ((St0)-max(D))^2 + ((Sp0)-min(D))^2;

if plot_flag, disp([abs(St0) max(D) abs(Sp0) min(D) res]);end
