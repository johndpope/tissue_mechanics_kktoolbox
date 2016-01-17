function X_o = invert_shp(X_o)


[xc yc zc] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);

for ix = 2:2:L_max      % only invert the sign of even L coefficients
    start = (ix-1 + 1)^2+1; 
    finish = (ix + 1).^2;
    vec = start:finish;
    xc(vec) = -xc(vec);
    yc(vec) = -yc(vec);
    zc(vec) = -zc(vec);
end
X_o = [xc(:)' yc(:)' zc(:)'];