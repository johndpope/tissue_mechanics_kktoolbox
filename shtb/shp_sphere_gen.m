function X_o = shp_sphere_gen(L_max,R)
% generates the set of coefficients that correspond to a sphere
% which means filling in the L = 1 coefficients
xclks = zeros(4,1);
yclks = zeros(4,1);
zclks = zeros(4,1);

xclks(4) = R/N_LK(1,1);
yclks(2) = R/N_LK(1,-1);
zclks(3) = R/N_LK(1,0);

X_o = [xclks(:)' yclks(:)' zclks(:)'];
X_o = tr(X_o, L_max);