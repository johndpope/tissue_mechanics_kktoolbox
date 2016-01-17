function ix = get_LK_index_P(L,K)
% assuming arrangement of indices L and K for associated Legendre funcions
% i.e. P00, P10, P11, P20 P21 ... etc., we get the linear
% index reflecting the position of the L,K pair
K = abs(K);
ix = uint8((L-1)^2/2+3*(L-1)/2 + 1) + K + 1;