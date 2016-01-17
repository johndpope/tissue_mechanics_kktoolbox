function [N, L, K] = sh4d_indices_gen(n,l)
% generate the indices in a consistent sequence for use with the 4D sh
% analysis
N = []; L = []; K = [];
[lvec, kvec] = shp_surface.indices_gen(1:(l + 1)^2);

%%% for the N = 0 case
    N = [N ones(1,length(lvec))*0];
    L = [L lvec(:)'];
    K = [K kvec(:)'];

for ix = 1:n
    %%% for cos
    N = [N ones(1,length(lvec))*ix];
    L = [L lvec(:)'];
    K = [K kvec(:)'];
    %%%uncomment for sin
    N = [N ones(1,length(lvec))*-ix];
    L = [L lvec(:)'];
    K = [K kvec(:)'];
end
N = N';L = L';K = K';
