function X_o = force_canonical2relaxed(X_o)
%%% only retains the major axes of the first order ellipsoid and the all
%%% coefficients with L > 1
%%% counterpart is force_canonical2relaxed

nc = round(length(X_o)/3);
xc = X_o(1:nc);
yc = X_o(nc+1:2*nc); 
zc = X_o(2*nc+1:3*nc);

xc = [0 0 0 xc(:)'];
yc = [0 yc(1) 0 0 yc(2:end)];
zc = [0 0 zc(1) 0 zc(2:end)];
X_o = [xc(:)' yc(:)' zc(:)'];