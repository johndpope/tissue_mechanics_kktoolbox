function Xtr = shp_truncation(X_o, fac)
% truncates the series by setting to zero coefficients below fac factor of the highest abs coeffients of
% the corresponding coordinate.
% Example: 
%       Xnew = shp_truncation(X_o, 0.01);
% Co;yright: Khaled Khairy 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xc yc zc] = get_xyz_clks(X_o);
xc(1) = 0;
yc(1) = 0;
zc(1) = 0;
xc(abs(xc)<(max(abs(xc))*fac))=0;
yc(abs(yc)<(max(abs(yc))*fac))=0;
zc(abs(zc)<(max(abs(zc))*fac))=0;

[xc2 yc2 zc2] = get_xyz_clks(X_o);
xc(1) = xc2(1);
yc(1) = yc2(1);
zc(1) = zc2(1);
Xtr = [xc(:)' yc(:)' zc(:)'];
disp(['Number of non-zero coefficients ' num2str(sum(Xtr~=0))]);
