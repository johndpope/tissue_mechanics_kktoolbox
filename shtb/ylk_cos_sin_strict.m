function [Y NLK P_LK] = ylk_cos_sin_strict(L,K,phi, theta)
% here we calclulate a real combination of the YLK's for  a specific value of L and K
% this function includes the Condon-Shortly phase factor introduced through
% Matlab\s lengendre function, and also for the case of negative K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if K ==0 && L ==0,
	Y = 1/2/sqrt(pi)* ones(size(theta));
else
	NLK = N_LK((L),(K));
	CS = (-1)^K;
    
    P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	P_LK = reshape(P_LK,size(theta,1),size(theta,2));   % Matlab already has the CS factor, so we leave P_LK as is.
    
    if K<0, 
        P_LK = P_LK * CS*factorial(L-abs(K))/factorial(L+abs(K)); % include CS factor when calculating PLK for negative K.
    end
    
	if K >= 0 
		Y = NLK * P_LK.*cos(K * phi);     %%% 
	else	% K is negative
		Y = NLK * P_LK.*sin(abs(K) * phi);
	end
end
