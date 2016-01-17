function Y = ylk_cos_sin(L,K,phi, theta)
% No Condon-Shortley phase factor version (the new version)
% here we calclulate a real combination of the YLK's for  specific values of L and K
% Since Matlab includes the Condon Shortley phase factor in its calculation
% of the associated Legendre polynomials, we are negating them for our
% calculation in compatibility with the functions that calcuate the
% derivatives of the Legendre polynomials recursively as in Duncan and
% Olson.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K ==0 && L ==0,
	Y = 1/2/sqrt(pi)* ones(size(theta));
else
	NLK = N_LK(L,K);%sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
	P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	P_LK = reshape(P_LK,size(theta,1),size(theta,2));
    CS = (-1)^K;
	if K >= 0 
		Y = CS*NLK * P_LK.*cos(K * phi);     %%% 
	else	% K < 0
		Y = CS*NLK * P_LK.*sin(abs(K) * phi);
	end
end
