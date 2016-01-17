function Coeff = DNLK_trunc_fixed(DNLK, n)
%% sorts the coefficients in ascending order and truncates

D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;

dim = length(D);
N = N(1:dim);L = L(1:dim);K = K(1:dim);

[ds, IX] = sort(D.^2, 'descend');


Coeff.D = ds(1:n);
Coeff.N = N(IX(1:n));
Coeff.L = L(IX(1:n));
Coeff.K = K(IX(1:n));