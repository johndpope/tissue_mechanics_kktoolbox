function X_o = axisymmetric2full(X_a)
%% turns the axisymmetric form X_a into full form. The shape is of course
%% still axisymmtric. It is assumed that X_a was generated with
%% full2axisymmetric.
L_list = 1:40;                   % let's neglect L = 0;
k0_ixs = L_list.^2 + L_list + 1;    % indices with k = 0;

max_k0 = k0_ixs(length(X_a) - 1);
L_max = (-1 + sqrt(1-4+4*max_k0))/2;

xclks = zeros((L_max + 1)^2,1);
xclks(4) = X_a(1);

yclks = zeros((L_max + 1)^2,1);
yclks(2) = X_a(1)*N_LK(1,1)/N_LK(1,-1);

zclks = zeros((L_max + 1)^2,1);
zclks(k0_ixs(1:length(X_a)-1)) = X_a(2:end);

X_o = [xclks(:)' yclks(:)' zclks(:)'];
