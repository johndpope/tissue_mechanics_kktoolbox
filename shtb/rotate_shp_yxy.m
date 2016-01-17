function S_rot = rotate_shp_yxy(S,ang, c)
%%% rotate the shapes described by the spherical harmonic coefficiencts S
%%% by the Euler angles a b g (radians) around "center"
%%% rotation conventions are y-x-y
c = c/sqrt(1/4/pi);
S_rot = [];
%% generate rotation matrix
a = ang(1);b = ang(2);g = ang(3);
Rg = rot_mx(g,2);
Rb = rot_mx(b,3);
Ra = rot_mx(a,2);
R = Ra*Rb*Rg;

%% loop over the shapes and rotate them around the center
for ix = 1:size(S,1),
    X_o = S(ix,:);
    [xc, yc, zc] = get_xyz_clks(X_o);

    xc(1) = xc(1)-c(1);
    yc(1) = yc(1)-c(2);
    zc(1) = zc(1)-c(3);

    C = [xc(:)';yc(:)';zc(:)'];
    cr = R*C;
    cr = cr';
    tx = [cr(:,1)];tx(1) = tx(1) + c(1);
    ty = [cr(:,2)];ty(1) = ty(1) + c(2);
    tz = [cr(:,3)];tz(1) = tz(1) + c(3);

    S_rot(ix,:) = [tx; ty; tz ];
end




