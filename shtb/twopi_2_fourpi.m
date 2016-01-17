function Xnew = twopi_2_fourpi(X_o)

nc = length(X_o)/3;yclks = X_o(nc+1:2*nc);zclks = X_o(2*nc+1:end);xclks = X_o(1:nc);
L_max = round(sqrt(length(xclks))-1);
counter = 1;
for L = 1:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
        NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
        NLK = N_LK(L,K);
        xclks(counter) = xclks(counter) * CS.*(NLK/NLK_old);
        yclks(counter) = yclks(counter) * CS.*(NLK/NLK_old);
        zclks(counter) = zclks(counter) * CS.*(NLK/NLK_old);
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];