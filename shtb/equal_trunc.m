function [X1tr, X2tr] = equal_trunc(X1,X2)
%%% produce shape vectors of equal Lmax (i.e. truncation)
Lmax1 = sqrt(length(X1)/3)-1;
Lmax2 = sqrt(length(X2)/3)-1;
L_max = max([Lmax1 Lmax2]);

clks = X1;nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
% if length(xclks)>(L_max+1)^2, xclks = xclks(1:(ud.L_max+1)^2);yclks = yclks(1:(ud.L_max+1)^2);zclks = zclks(1:(ud.L_max+1)^2);end
if length(xclks)<(L_max+1)^2, xclks((L_max+1)^2) = 0;yclks((L_max+1)^2) = 0;zclks((L_max+1)^2) = 0;end
X1tr = [xclks(:)' yclks(:)' zclks(:)'];

clks = X2;nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
% if length(xclks)>(L_max+1)^2, xclks = xclks(1:(ud.L_max+1)^2);yclks = yclks(1:(ud.L_max+1)^2);zclks = zclks(1:(ud.L_max+1)^2);end
if length(xclks)<(L_max+1)^2, xclks((L_max+1)^2) = 0;yclks((L_max+1)^2) = 0;zclks((L_max+1)^2) = 0;end
X2tr = [xclks(:)' yclks(:)' zclks(:)'];