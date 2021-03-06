function Y_TP = precalc_ylk_cos_sin_dthetaphi(phi, theta, L_max)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
global P_LK_T
gdimp = size(phi,2);
gdimt = size(theta,1);

if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:);
    theta = theta(:);
    Y_TP = zeros(gdimt,1, (L_max+1)^2);
else
    Y_TP 	= zeros(gdimt, gdimp, (L_max+1)^2);
end



counter 	= 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        NLK = N_LK(L,K);%NLK = sqrt((2 * L + 1)/(4*pi) * factorial(L - abs(K))/factorial(L + abs(K)));
        if K == 0,
            Y_TP(:,:,counter) = zeros(size(Y_TP(:,:,counter)));    
        elseif K > 0,
            Y_TP(:,:,counter) =  NLK * P_LK_T(:,:, counter).* (-K).*sin(K * phi);%
        elseif K < 0,
            Y_TP(:,:,counter) =  NLK * P_LK_T(:,:, counter).* abs(K).* cos(abs(K) * phi);%
        end
     end
end



