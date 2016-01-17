function Y_P = ylk_cos_sin_dphi_bosh(p, t, L_max,P)
%%% calculate the derivative of the normalized associated Legendre functions
%%% up to L = L_max
%%% needs P as input (obtained as output from ylk_cos_sin_bosh

gdimp = length(p);
gdimt = length(t);
if (min(size(p))==1 && max(size(t))==1),
     p = p(:);
end
    Y_P = zeros(gdimp, gdimt, (L_max+1)^2); 

counter = 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        if K == 0,
            Y_P(:,:,counter) = zeros(size(Y_P(:,:,counter)));
        elseif K > 0,
            Y_P(:,:,counter) =  P(:,:,get_LK_index_P(L,K)).* (-K).*sin(K.*p);%
        elseif K < 0,
            Y_P(:,:,counter) =  P(:,:,get_LK_index_P(L,K)).*abs(K).*cos(abs(K).*p);%
        end
     end
end




