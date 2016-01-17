function X_o = fix_xc20(X_o)


[xc yc zc] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);
last = 0;
for ix = 0:L_max
    start = last+1; finish = (ix + 1).^2;
    vec = start:finish;
    
    last = finish;
    if kk_iseven(ix),    % i.e. if ix is odd
        xc(vec) = -xc(vec);
        yc(vec) = -yc(vec);
        zc(vec) = -zc(vec);
    end
end
X_o = [xc(:)' yc(:)' zc(:)'];


