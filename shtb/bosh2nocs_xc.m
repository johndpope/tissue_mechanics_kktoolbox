function xclks = bosh2nocs_xc(xclks)

L_max = round(sqrt(length(xclks))-1);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        NLK_old = sh_basis.N_LK_nocs(L,K);
        NLK = sh_basis.N_LK_bosh(L,K);
        fac = NLK_old/NLK;
        xclks(counter) = xclks(counter) * 1/fac;
    end
end
