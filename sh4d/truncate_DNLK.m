function DNLK = truncate_DNLK(DNLK, thresh)
%% truncate a DNLK coefficient set according to the amplitude of the D
%% coefficients
%% thresh = factor of maximum coefficient squared

D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;
tau = DNLK.tau;
disp('Old basis set dimension :');disp(length(D));

D = D(D.^2>(max(D.^2)*thresh));
N = N(D.^2>(max(D.^2)*thresh));
L = L(D.^2>(max(D.^2)*thresh));
K = K(D.^2>(max(D.^2)*thresh));

DNLK.D = D;
DNLK.N = N;
DNLK.L = L;
DNLK.K = K;


disp('New basis set dimension :');disp(length(D));