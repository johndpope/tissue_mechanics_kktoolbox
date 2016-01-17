function X_1 = restore_configuration_from_r_inv(X_r,J)
%% rotate and invert the shape described by the spherical harmonic coefficiencts X_r
%% after it has been subjected to rotation in object space (T_2) and possible
%% inversion operations
%% J is the vector of operations (three Euler angles and inversion) as
%% obtained from r_inv, that will be used to bring back the shape configuration to look like the original X_o and
%% to be equal to the parametric-canonical shape X_1;

[xc yc zc] = get_xyz_clks(X_r);

xc = xc*J(4);yc = yc*J(5);zc = zc*J(6);
X_r = [xc(:)' yc(:)' zc(:)'];

ang = J(1:3);
X_1 = rotate_shp(X_r, -(ang(3:-1:1)));

