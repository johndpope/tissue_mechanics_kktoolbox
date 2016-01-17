function xclks = trunc_sh(clks, L_max)


    nc = round(length(clks));
    xclks = clks(1:nc);
    trunc = (L_max+1)^2;
    lmax_in = sqrt(length(clks))-1;
    if lmax_in<L_max,
        xclks(trunc) = 0;
    else
        xclks = xclks(1:trunc);
    end
    
