function [Y_LK_theta, P_LK_T] = precalc_ylk_cos_sin_dtheta_ng(phi, theta, L_max, P_LK)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max

gdimp = size(phi,2);
gdimt = size(theta,1);

if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:);
    theta = theta(:);
    Y_LK_theta = zeros(gdimt,1, (L_max+1)^2);
    P_LK_T = Y_LK_theta;
else
Y_LK_theta 	= zeros(gdimt, gdimp, (L_max+1)^2); 
P_LK_T 		= zeros(gdimt, gdimp, (L_max+1)^2); 
end




counter 	= 0;%
for L = 0:L_max
    [p plk_theta] = plktheta(L,cos(theta), theta, 0, P_LK);   
    for K = -L:L
        counter = counter + 1;
        NLK = N_LK(L,K);%NLK = sqrt((2 * L + 1)/(4*pi) * factorial(L - abs(K))/factorial(L + abs(K)));
        if K >= 0,
            Y_LK_theta(:,:,counter) = NLK * plk_theta(:,:, K+1).* cos(K * phi);%
            P_LK_T(:,:,counter) 	=  plk_theta(:,:, K+1);%
        elseif K < 0,
            pK = abs(K);
            Y_LK_theta(:,:,counter) = NLK * plk_theta(:,:,abs(K)+1).* sin(abs(K) * phi);%
            P_LK_T(:,:,counter) 	= plk_theta(:,:, abs(K)+1);%
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [p plk_theta] = plktheta(L,ct, theta, K, P_LK, plk_theta)
%%% Calculate the derivative of the associated legendre functions using the
%%% recursion relations from Duncan and Olson '93

if K == L,
    if K == 0 ,
%         plk_theta = zeros(length(theta), length(theta), L+1);
        plk_theta = zeros(size(P_LK,1), size(P_LK,2), L+1);
    end
%    p = (-1)^K *L * factorial(2*L)/2/factorial(L) * (sin(theta)/2).^(L-1) .* ct;
    p = L .* factorial(2*L)/2/factorial(L) .* (sin(theta)/2).^(L-1) .* ct;

    plk_theta(:,:,K+1) = p;
elseif K > L,
    p = zeros(size(P_LK,1), size(P_LK,2));
elseif K < L,
    if K == 0 ,
%         plk_theta = zeros(length(theta), length(theta), L+1);
        plk_theta = zeros(size(P_LK,1), size(P_LK,2), L+1);
    end
    [plktp1 plk_theta] = plktheta(L,ct,theta,K+1, P_LK, plk_theta);
    [plktp2 plk_theta]  = plktheta(L,ct,theta, K+2, P_LK, plk_theta);
%     plkp1  = (-1)^K*P_LK(:,:,(L+1)^2 - (L-(K+1)));
     plkp1  = P_LK(:,:,(L+1)^2 - (L-(K+1)));
%    plkp1  = P_LK(:,:,K+2);
    
    p = 1/(L+K+1)/(L-K) * ...
        ( 2*(K+1).*cot(theta).* plktp1 - 2*(K+1).*csc(theta).^2.*plkp1-plktp2);
    plk_theta(:,:,K+1) = p;    
end




