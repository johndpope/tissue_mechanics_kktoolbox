function Y = ylk_bosh(L,K,phi, theta)
% here we calclulate a real combination of the YLK's for  a specific value of L and K
% this function negates the Condon-Shortly phase factor introduced through
% Matlab's lengendre function, and uses the Bosh 2000 normalization
%%% Author: Dr. Khaled Khairy (khaledkhairy@yahoo.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if K ==0 && L ==0,
	Y = 1/2/sqrt(pi)* ones(size(theta));
else
	NLK = N_LK((L),(K));
	P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	P_LK = reshape(P_LK,size(theta,1),size(theta,2));
    CS = (-1)^K;
	if K >= 0 
		Y = CS*NLK_bosh * P_LK.*cos(K * phi);     %%% 
	else	% K < 0
		Y = CS*NLK_bosh * P_LK.*sin(abs(K) * phi);
	end
end
