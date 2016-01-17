function S_rot = rotate_x_shps(S,ang, c)
%%% rotate the shapes described by the spherical harmonic coefficiencts S
%%% around the x-axis around the center c
c = c/sqrt(1/4/pi);
S_rot = [];
%% generate rotation matrix
Rx = rot_mx(ang,3);
R = Rx;
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




