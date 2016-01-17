function X_a = full2axisymmetric(X_o)
%% returns a shortened X_o that represents an axisymmetric shape with the
%% axis of symmetry z and centered at zero. The vector's first value X_a(1)
%% is xclk(4), and X_a(2) corresponds to yclk(2). There is no L0K0 entries.

L_max = get_L_max(X_o);
[xclks , ~, zclks] = get_xyz_clks(X_o);
xy = xclks(4);
L_list = 1:L_max;                   % let's neglect L = 0;
k0_ixs = L_list.^2 + L_list + 1;    % indices with k = 0;
%tempy = tempx*N_LK(1,1)/N_LK(1,-1);
X_a = [xy zclks(k0_ixs)'];
