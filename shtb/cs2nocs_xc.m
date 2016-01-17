function xclks = cs2nocs_xc(xclks)
L_max = round(sqrt(length(xclks))-1);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
        NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
        NLK = sh_basis.N_LK_nocs(L,K);
        xclks(counter) = xclks(counter) * CS*NLK_old/NLK;
    end
end
