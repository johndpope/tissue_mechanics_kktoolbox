function X_o = set_axisymmetric(X_o)
L_max = sqrt(round(length(X_o)/3))-1;
[xclks yclks zclks] = get_xyz_clks(X_o);
tempx = xclks(4);
tempz = zclks(3);
L_list = 1:L_max;                   % let's neglect L = 0;
k0_ixs = L_list.^2 + L_list + 1;          % indices with k = 0;
xclks(:) = 0;xclks(4) = tempx;%temp = xclks;xclks(k0_ixs) = temp(k0_ixs);
yclks(:) = 0;yclks(2) = tempx*N_LK(1,1)/N_LK(1,-1);%temp = yclks;yclks(k0_ixs) = temp(k0_ixs);
temp = zclks;zclks(:) = 0;zclks(k0_ixs) = temp(k0_ixs); zclks(3) = tempz;

X_o = [xclks(:)' yclks(:)' zclks(:)'];
