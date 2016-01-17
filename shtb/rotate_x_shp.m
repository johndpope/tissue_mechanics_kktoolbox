function X_rot = rotate_x_shp(X_o,ang)
%% rotate the shape described by the spherical harmonic coefficiencts X_o
%% around x axis by angle (radians). Returns the rotated coefficients

[xc, yc, zc] = get_xyz_clks(X_o);
%% rotation conventions are y-z-y
x = xc(1);xc(1) = 0;
y = yc(1);yc(1) = 0;
z = zc(1);zc(1) = 0;

Ra = rot_mx(ang,3);
C = [xc(:)';yc(:)';zc(:)'];
cr = [Ra*C];
cr = cr';
tx = [cr(:,1)];tx(1) = x;
ty = [cr(:,2)];ty(1) = y;
tz = [cr(:,3)]; tz(1) = z;
X_rot = [tx; ty; tz ]; 



