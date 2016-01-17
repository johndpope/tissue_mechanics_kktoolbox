function [Y_LK, P_LK NLK_vec]= precalc_ylk_cos_sin(phi, theta, L_max)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max... note this version is the new one without the CS phase
%% factor. and the new NLK normalization which makes it compatible with the
%% functions that calculate the derivatives later on.
global plk
%plk = [];
%%%%% The following has been tested and works!
gdimp = size(phi,2);
gdimt = size(theta,1);
NLK_vec = [];
if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:);
    theta = theta(:);
    Y_LK = zeros(gdimt,1, (L_max+1)^2);
    P_LK = Y_LK;
else
    Y_LK = zeros(gdimt, gdimp, (L_max+1)^2);
    P_LK = zeros(gdimt, gdimp, (L_max+1)^2);
end
counter = 0;%
for L = 0:L_max
    plk_mat = legendre(L,cos(theta(:)'));
    for K = -L:L
        counter = counter + 1;
        plk = reshape(plk_mat(abs(K)+1,:),size(theta,1),size(theta,2));
        NLK = N_LK(L,K);%sqrt((2 * L + 1)/(4*pi) * factorial(L - (K))/factorial(L + (K)));
        NLK_vec(counter) = NLK;
        %%%fac = (-1)^(K).*factorial(L-abs(K))./factorial(L+abs(K));
        CS = (-1)^K;
        if K ==0 && L ==0
            Y_LK(:,:,counter) =  NLK.*plk;
            P_LK(:,:,counter) =  plk;
        elseif K >=0,
            P_LK(:,:,counter) = CS*(plk);
            Y_LK(:,:,counter) = CS*NLK.* plk.*cos(K.*phi);
        elseif K<0,
            Y_LK(:,:,counter) = CS*NLK*plk.* sin(abs(K).*phi);
            %             Y_LK(:,:,counter) = fac*NLK.* squeeze(plk(abs(K)+1,:,:)).* sin(abs(K).*phi);
            P_LK(:,:,counter) = CS*plk;
            %            P_LK(:,:,counter) =  factorial(L-pK)./factorial(L+pK).*squeeze(plk(abs(K)+1,:,:));
        end
    end
end
