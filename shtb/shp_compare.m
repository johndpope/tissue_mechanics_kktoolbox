function shp_compare(X1, X2)
%%% For debugging purposes let us look at the coefficients of 2 shapes in
%%% some detail.
%%% Input X1 and X2 need to be rotation invariant as if resulting from
%%% r_inv. They do not need to be scale invariant (this is done here).

L_max = 1;
X1 = s_inv(tr(X1,L_max));
X2 = s_inv(tr(X2,L_max));
[xc yc zc] = get_xyz_clks(X1);
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD1 = sqrt(diag(D));


[xc1 yc1 zc1] = get_xyz_clks(X2);
Ar = [xc1(2:4) yc1(2:4) zc1(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD2 = sqrt(diag(D));
disp([xc(:) xc1(:) yc(:) yc1(:) zc(:) zc1(:)]);

% % load clk_str;
% % for ix = 1:min([length(xc) length(clk_str)])
% % disp([clk_str{ix} '    ' num2str([xc(ix)-xc1(ix) yc(ix)-yc1(ix) zc(ix)-zc1(ix)])]);
% % end
disp('Diagonal');
disp(sqrtdD1);disp('---');disp(sqrtdD2);disp('----');
D = [xc(:)-xc1(:) yc(:)-yc1(:) zc(:)-zc1(:)].^2;

% % figure;
% % subplot(3,1,1);plot(D(1,:));
% % subplot(3,1,2);plot(D(2,:));
% % subplot(3,1,3);plot(D(3,:));