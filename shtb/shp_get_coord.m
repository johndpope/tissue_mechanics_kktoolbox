function [X,C]=shp_get_coord(X_o, gdim)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = round(sqrt(length(xclks))-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];