function Y = ylk_cos_sin_old(L,K,phi, theta)
% Old version. This function is not compatible with the expressions for
% calculating the derivatives of the Legendre Polynomials, since the
% Condon-Shortley phase factor which matlab includes in its evaluation of
% the associated Legendre functions was not negated. Also note that the
% normalization NLK here is different from the new version.
% here we calclulate a real combination of the YLK's for  a specific value of L and K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K ==0 & L ==0,
	Y = 1/2/sqrt(pi)* ones(size(theta));
else
	NLK = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
	P_LK = legendre(L,cos(theta));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	
	if K >= 0 
		Y = NLK * P_LK.*cos(K * phi);     %%% 
	else	% K < 0
		Y = NLK * P_LK.*sin(abs(K) * phi);
	end
end
