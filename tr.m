function S_tr = tr(S, L_max)

if size(S,2) ==1,
    clks = S(:)';% check whether S has a dimension which is one. i.e. we have only one shape
    nc = round(length(clks)/3);
    xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
    trunc = (L_max+1)^2;
    lmax_in = sqrt(length(clks)/3)-1;
    if lmax_in<L_max,
        xclks(trunc) = 0;
        yclks(trunc) = 0;
        zclks(trunc) = 0;
    else
        xclks = xclks(1:trunc);
        yclks = yclks(1:trunc);
        zclks = zclks(1:trunc);
    end
    S_tr = [xclks(:)' yclks(:)' zclks(:)']';
else
    S_tr = zeros(size(S,1), 3*(L_max+1)^2);
    for ix = 1:size(S,1),       % loop over the shapes
        S_tr(ix,:) = tr(S(ix,:)', L_max);
    end
end