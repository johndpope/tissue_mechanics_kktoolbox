function [all_X, all_F] = get_all_X_from_T(T,gdim)
%% returns the cell array of coordinates after reconstucting the object
%% surface from the SHP
L_max = round(sqrt(size(T,2)/3)-1);

P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

Y_LK = get_basis(t',p',gdim,L_max);
all_X = {};
all_F = {};

for ix = 1:size(T,1)        % loop over objects
    P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
    [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
    X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
    [xclks yclks zclks] = get_xyz_clks(T(ix,:));
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
    all_X{ix} = X;
    all_F{ix} = C;
end