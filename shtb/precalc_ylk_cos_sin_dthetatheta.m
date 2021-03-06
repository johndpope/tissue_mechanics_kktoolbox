function [Y_TT, P_TT] = precalc_ylk_cos_sin_dthetatheta(phi, theta, L_max)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
global ptt
gdimp = size(phi,2);
gdimt = size(theta,1);
[tt wtt]                  = gaussquad(gdimt, 0, pi);

if (min(size(phi))==1 && min(size(theta))==1),
    phi = phi(:);
    theta = theta(:);
    Y_TT = zeros(gdimp,1, (L_max+1)^2);
    P_TT = Y_TT;
else

    Y_TT = zeros(gdimt, gdimp, (L_max+1)^2);
    P_TT = zeros(gdimt, gdimp, (L_max+1)^2);
end

 
counter = 0;%
for L = 0:L_max
        p = ptt_gen(L,cos(theta), theta, 0);   
    for K = -L:L
        counter = counter + 1;
        NLK = N_LK(L,K);%NLK = sqrt((2 * L + 1)/(4*pi) * factorial(L - abs(K))/factorial(L + abs(K)));
        if K >= 0,
%              if K==0,
%                ytt = ptt(:,:,K+1);
%                 ytt(theta==(tt(1)))= ytt(theta==tt(4));
%                 ytt(theta==(tt(2)))= ytt(theta==tt(4));
%                 ytt(theta==(tt(3)))= ytt(theta==tt(4));
%                 
%                 ytt(theta==(tt(end))) = ytt(theta==tt(end-3));
%                 ytt(theta==(tt(end-1))) = ytt(theta==tt(end-3));
%                 ytt(theta==(tt(end-2))) = ytt(theta==tt(end-3));
%                 ptt(:,:,K+1) = ytt;
%             end
            Y_TT(:,:,counter) =  NLK * ptt(:,:, K+1).* cos(K.*phi);%
% 			P_TT(:,:,counter) =  ptt(:,:, K+1);%
        elseif K < 0,
%             pK = abs(K);
            Y_TT(:,:,counter) =  NLK.* ptt(:,:,abs(K)+1).* sin(abs(K).* phi);%
%    			P_TT(:,:,counter) =  (-1)^pK.* factorial(L-pK)./factorial(L+pK).*ptt(:,:,abs(K)+1);%
        end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function p = ptt_gen(L,ct, theta, K)
%%% Calculate the derivative of the associated legendre functions using the
%%% recursion relations from Duncan and Olson '93
global P_LK P_LK_T ptt
if K == L,
    if K == 0 ,
        ptt = zeros(size(P_LK,1), size(P_LK,2), L+1);
    end
%     p = (-1)^K * -L * factorial(2*L)/4/factorial(L) * (sin(theta)/2).^(L-2) .* (1-L * ct.^2);

    p =  -L * factorial(2*L)/4/factorial(L) * (sin(theta)/2).^(L-2) .* (1-L * ct.^2);

    ptt(:,:,K+1) = p;
elseif K > L,
    p = zeros(size(P_LK,1), size(P_LK,2));
elseif K < L,
    if K == 0 ,
        ptt = zeros(size(P_LK,1), size(P_LK,2), L+1);
    end
    pttp1 = ptt_gen(L,ct,theta,K+1);
    pttp2 = ptt_gen(L,ct,theta, K+2);
    plktp1 = P_LK_T(:,:,(L+1)^2 - (L-(K+1)));
    plkp1 = P_LK(:,:,(L+1)^2 - (L-(K+1)));
    
% %     fac = (L+K+1)/(L-K);
% %     p = ...
% %         1/fac*(...
% %         -2.*(2.*K+1).*cot(theta).*(-1-cot(theta).^2).*plkp1+...
% %             2.*(2.*K+1).*(-1-cot(theta).^2).*plktp1+...
% %             (2.*K+1).*cot(theta).*pttp1-...
% %             pttp2...
% %         );
    
    
    p = 1/(L+K+1)/(L-K) * ((2*(K+1)-1).*...
        (cot(theta).*pttp1 ...
        - 2 .*csc(theta).^2 .* plktp1 ...
		+ 2 .* cot(theta).*csc(theta).^2 .* plkp1)...
        - pttp2);

%     p = 1/(L+K+1)/(L-K) * (2 * (K+1) .*(cot(theta).*pttp1 ...
%         - 1 .*csc(theta).^2 .* plktp1 ...
% 		+ 2 .* cot(theta).*csc(theta).^2 .* plkp1)...
%         - pttp2);
    ptt(:,:,K+1) = p;
end




