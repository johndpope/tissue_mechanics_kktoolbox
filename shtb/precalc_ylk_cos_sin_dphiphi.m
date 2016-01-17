function Ypp = precalc_ylk_cos_sin_dphiphi(phi, theta, L_max)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
global P_LK
gdimp = size(phi,2);
gdimt = size(theta,1);
if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:);
    theta = theta(:);
    Ypp = zeros(gdimt,1, (L_max+1)^2);

else
    Ypp = zeros(gdimt, gdimp, (L_max+1)^2); 
end


counter = 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        NLK = N_LK(L,K);
        if K == 0,
            Ypp(:,:,counter) = zeros(size(Ypp(:,:,counter)));
        elseif K > 0,
            Ypp(:,:,counter) = NLK.* P_LK(:,:,counter).*-(K^2).*cos(K.* phi);%
        elseif K < 0,
            Ypp(:,:,counter) = NLK.* P_LK(:,:,counter).*-(K^2).*sin(abs(K).* phi);%
        end
     end
end




