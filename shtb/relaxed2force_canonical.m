function X_f = relaxed2force_canonical(X_o)
%%% only retains the major axes of the first order ellipsoid and the all
%%% coefficients with L > 1
%%% counterpart is force_canonical2relaxed

[xc yc zc] = get_xyz_clks(X_o);
xc([1 2 3]) = [];
yc([1 3 4]) = [];
zc([1 2 4]) = [];
X_f = [xc(:)' yc(:)' zc(:)'];