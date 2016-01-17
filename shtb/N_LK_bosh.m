function NLK = N_LK_bosh(L,K)
%%% calculates the normalization factor compatible with the Bosch expressions
%%% Author: Dr. Khaled Khairy (khaledkhairy@yahoo.com)
K = abs(K);

if abs(K)>L, 
	NLK = 0;
else
	NLK = sqrt((2-isequal(K,0))*(2*L+1)*factorial(L-K)/factorial(L+K));
end
