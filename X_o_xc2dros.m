function d = X_o_xc2dros(X_o,xc)
%% convert shape encoded in X_o and field ecoded in xc into a dros_embryo object
L_max = shp_surface.get_L_max(X_o);
gdim = 180;
basis = sh_basis(L_max,gdim);

d = dros_embryo(L_max,basis);
d.X_o = X_o;
s = sh_surface(L_max,basis);
s.xc = xc;
d.sf{1} = {'untitled', s};
