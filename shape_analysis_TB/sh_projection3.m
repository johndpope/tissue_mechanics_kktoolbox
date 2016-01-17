function [X_o] = sh_projection3(L_max, X, t, p)
%%% The expansion of the three functions x(t,p), y(t,p) and z(t,p) on the sphere.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, K ] = indices_gen(1:(L_max + 1)^2);
M = length(L); %% number of functions in expansion
N = size(X,1); %% number of data points
A  = zeros(N, M, 'single');
for S = 1:length(L),  A(:,S) = ylk_cos_sin(L(S),K(S),p(:)',t(:)')'; end% prepare basis functions
[U, S, V] = svd(A, 'econ');
warning off;invS = 1./(S);invS(invS==inf) = 0;warning on
xclks = (V*invS) * (U'*X(:,1));
yclks = (V*invS) * (U'*X(:,2));
zclks = (V*invS) * (U'*X(:,3));
X_o = [yclks(:)' zclks(:)' xclks(:)'];




