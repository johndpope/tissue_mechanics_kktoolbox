function X_s = s_inv(X_o)
%%% Generate scale invariance based on the a normalization with respect to
%%% the longest axis of the first order ellipsoid. Translational invariance
%%% is also done in this function.

[xc yc zc] = get_xyz_clks(X_o);
xc(1) = 0;yc(1) = 0;zc(1) = 0;
X_s = [xc(:)' yc(:)' zc(:)']./abs(zc(3));