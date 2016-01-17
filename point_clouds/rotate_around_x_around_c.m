function X_r = rotate_around_x_around_c(X_o,a,c)
%%% returns the coefficients after rotation about a center y,z defined
%%% in C (a 3-vector) by angle a in degrees


[xco, yco, zco] = get_xyz_clks(X_o);
Ra = rot_mx(deg2rad(a),3);
% translate such that the center of mass is at the rotation axis
xc = xco;yc = yco;zc = zco;
%xc1 = xc;xc1(1) = xc(1) - c(1);
fac = sqrt(1/4/pi);
yc(1) =  yc(1) - c(2)/fac;
zc1(1) = zc(1) - c(3)/fac;
% apply the rotation
C = [xc(:)';yc(:)';zc(:)'];
cr = [Ra*C];
cr = cr';


tx = [cr(:,1)];tx(1)  = xco(1);
ty = [cr(:,2)];ty(1)  = ty(1)+c(2)/fac;%yc(1);
tz = [cr(:,3)];tz(1) = tz(1)+c(3)/fac;%zc(1);
X_r = [tx; ty; tz ];