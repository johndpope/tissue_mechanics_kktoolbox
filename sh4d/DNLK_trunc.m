function DNLK = DNLK_trunc(DNLK, s)
%% tuncate the series by using only the first s number of basis vectors

D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;

DNLK.D = D(1:s);
DNLK.N = N(1:s);
DNLK.L = L(1:s);
DNLK.K = K(1:s);
