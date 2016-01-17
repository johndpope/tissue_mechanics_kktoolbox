function DNLK = DNLK_trunc_K(DNLK, s)
%% tuncate the series by using only the K = s coefficients

DNLK.old = DNLK;
temp = DNLK;
DNLK.D = [];
DNLK.N = [];
DNLK.L = [];
DNLK.K = [];

for ix = 1:length(s),
    D = temp.D;
    N = temp.N;
    L = temp.L;
    K = temp.K;

    indx = find(K==s(ix));
    D = D(indx);
    N = N(indx);
    L = L(indx);
    K = K(indx);

    DNLK.D = [DNLK.D;D];
    DNLK.N = [DNLK.N;N];
    DNLK.L = [DNLK.L;L];
    DNLK.K = [DNLK.K;K];
end