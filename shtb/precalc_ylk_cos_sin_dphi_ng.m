function Y_LK_phi = precalc_ylk_cos_sin_dphi(phi, theta, L_max, P_LK)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
gdimp = size(phi,2);
gdimt = size(theta,1);

if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:); theta = theta(:);
    Y_LK_phi = zeros(gdimt,1, (L_max+1)^2);
else
    Y_LK_phi = zeros(gdimt, gdimp, (L_max+1)^2); 
end


counter = 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        NLK = N_LK(L,K);%sqrt((2 * L + 1)/(4*pi) * factorial(L - abs(K))/factorial(L + abs(K)));
        if K == 0,
            Y_LK_phi(:,:,counter) = zeros(size(Y_LK_phi(:,:,counter)));
        elseif K > 0,
            Y_LK_phi(:,:,counter) =  NLK .* P_LK(:,:,counter).* (-K).*sin(K.*phi);%
        elseif K < 0,
            Y_LK_phi(:,:,counter) =  NLK .* P_LK(:,:,counter).*abs(K).*cos(abs(K).*phi);%
        end
     end
end




