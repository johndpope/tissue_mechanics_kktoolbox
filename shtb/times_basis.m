function X_b = times_basis(X_o)
X_b = X_o(:);
% %%% factor out the basis normalization before calculating the rotation
% [L, K] = indices_gen(1:(get_L_max(X_o) + 1)^2);
% N = [];
% for S = 1:length(L),N(S) = N_LK_bosh(L(S),K(S)); end% usese nocs version
% X_b = X_o(:).*[N(:)' N(:)' N(:)']';
% %%%
% 
