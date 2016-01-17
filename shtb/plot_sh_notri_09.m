function [A, V,h,T, dA, E, H, wp, wt, p, t]  = plot_sh_notri_08(X_o, gdim)
% Calculates the shape properties of a shape defined using the spherical harmonics parameterization
% shape properties include:
%       Area, bending energy, volume, total curvature and Gaussian
%       curvature
% Author: Khaled Khairy 2010
%%%%%% XXX DOES NOT WORK YET %%%%%%%
verbose = 0;plotflag = 1;
if nargin == 1,
    gdim = 60;
end
[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);

[t wt]                  = gaussquad(gdim, 0, pi );
[p wp]                  = gaussquad(gdim,0, 2*pi);
[p t]                   = meshgrid(p,t);
[wp wt]                 = meshgrid(wp, wt);

[Y P]                   = ylk_cos_sin_bosh(p, t, L_max);
Y_P                     = ylk_cos_sin_dphi_bosh(p, t, L_max, P);
[Y_T P_T]               = ylk_cos_sin_dtheta_bosh(p, t, L_max, P);
Y_PP 					= ylk_cos_sin_dphiphi_bosh(p, t, L_max, P);
Y_TT 					= ylk_cos_sin_dthetatheta_bosh(p, t, L_max,P_T);
Y_TP 					= ylk_cos_sin_dthetaphi_bosh(p, t, L_max, P_T);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdimp = length(p);gdimt = length(t);
c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
x = sum(c.*Y,3);
xp =(sum(c.*Y_P,3));
xt = (sum(c.*Y_T,3));
xpp = (sum(c.*Y_PP,3));
xtt = (sum(c.*Y_TT,3));
xtp = (sum(c.*Y_TP,3));

c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
y = sum(c.*Y,3);
yp =(sum(c.*Y_P,3));
yt = (sum(c.*Y_T,3));
ypp = (sum(c.*Y_PP,3));
ytt = (sum(c.*Y_TT,3));
ytp = (sum(c.*Y_TP,3));

c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
z = sum(c.*Y,3);
zp =(sum(c.*Y_P,3));
zt = (sum(c.*Y_T,3));
zpp = (sum(c.*Y_PP,3));
ztt = (sum(c.*Y_TT,3));
ztp = (sum(c.*Y_TP,3));

%%%%%%%%%%%%%%%%%%%%%%%% Calculate the geometrically important quantities
wp = wp(:);wt = wt(:);w = wp.*wt;

% % %   xt(t>pi/2) = -xt(t>pi/2);
% % %   yt(t>pi/2) = -yt(t>pi/2);
% % %   zt(t>pi/2) = -zt(t>pi/2);
% % % % % 
% % %  xtp(t>pi/2) = -xtp(t>pi/2);
% % %  ytp(t>pi/2) = -ytp(t>pi/2);
% % %  ztp(t>pi/2) = -ztp(t>pi/2);

X  =[x(:) y(:) z(:)];
Xt = [xt(:) yt(:) zt(:)];
Xp = [xp(:) yp(:) zp(:)];
Xpp = [xpp(:) ypp(:) zpp(:)];
Xtp = [xtp(:) ytp(:) ztp(:)];
Xtt = [xtt(:) ytt(:) ztt(:)];




E = dot(Xt,Xt,2);
E_test = (cos(t(:))).^2 .*(16*(sin(p(:)).^2) + (cos(p(:))).^2) + (sin(t(:))).^2;

F = dot(Xt,Xp,2);
G = dot(Xp,Xp,2);
SS = (cross(Xt,Xp,2));SSn = sqrt(E.*G-F.*F);
n = SS./SSn(:,ones(1,3));
L = dot(Xtt,n,2);
M = dot(Xtp,n,2);
N = dot(Xpp,n,2);

%% calculation of geometric properties
V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
A = sum(w.*SSn);
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));

%% %
%H = (( (xt).*xt+ (yt).*yt+ (zt).*zt).*(- (xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))+( (xp).*xp+ (yp).*yp+ (zp).*zp).*(- (xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp))))-(2.* (xt).*xp+2.* (yt).*yp+2.* (zt).*zp).*(1./2.* (xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))))./(2.*( (xt).*xt+ (yt).*yt+ (zt).*zt).*( (xp).*xp+ (yp).*yp+ (zp).*zp)-2.*( (xt).*xp+ (yt).*yp+ (zt).*zp).^2);H = -H;

h = 1./A.*sum(H(:).*w.*SSn);
K = ((-conj(xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).*(-conj(xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))-(1./2.*conj(xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).^2)./((conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp)-(conj(xt).*xp+conj(yt).*yp+conj(zt).*zp).^2);
T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
dA= sum(sum(w.*2.*H(:).*SSn));         % assuming unit distance D (since D goes away in some expressions below)
% Curvedness = (2.*H.^2-K).^(0.5);		% curvedness
%fac = (H.^2-K).^(1/2);k_min = H-fac;k_max = H+fac;
energy = 1/2*sum((2.*H(:)).^2.*w.*SSn)./8/pi;   % reduced bending energy

if plotflag,
   %%
    dfig(1);clf;
    C = reshape(H,size(x));
    minH = -10 ;  C(C<minH) = minH;
    maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90)
     xlabel('x');ylabel('y');zlabel('z');colorbar
  ix = 1;
    dfig(2);clf;
    C = reshape(E,size(x));
%     minH = -10 ;  C(C<minH) = minH;
%     maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90)
    xlabel('x');ylabel('y');zlabel('z');title('E');
    
    dfig(3);clf;
    C = reshape(G,size(x));
%     minH = -10;  C(C<minH) = minH;
%     maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('G');
    
    dfig(4);clf;
    C = reshape(F,size(x));
%     minH = -10 ;  C(C<minH) = minH;
%     maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('F');

    dfig(5);clf;
    C = reshape(L,size(x));
    minH = -10 ;  C(C<minH) = minH;
    maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('L');
    
    dfig(6);clf;
    C = reshape(M,size(x));
    minH = -10 ;  C(C<minH) = minH;
    maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('M');
    
    dfig(7);clf;
    C = reshape(N,size(x));
    minH = -10 ;  C(C<minH) = minH;
    maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('N');

    dfig(8);clf;
    C = reshape(t,size(x));
    minH = -10 ;  C(C<minH) = minH;
    maxH = 10 ;  C(C>maxH) = maxH;
    surf(x,y,z,double(C), 'EdgeColor', 'none');axis equal;view(0,-90) 
    xlabel('x');ylabel('y');zlabel('z');title('t');
    
    dfig(11);clf;
    C = reshape(Xt(:,1),size(x));
    surf(p/pi, t/pi,double(C), 'EdgeColor', 'none');axis equal;view(0,-90);
    xlabel('phi');ylabel('theta');
    
    dfig(12);clf;
    C = reshape(Xt(:,2),size(x));
    surf(p/pi, t/pi,double(C), 'EdgeColor', 'none');axis equal;view(0,-90);
    xlabel('phi');ylabel('theta');
    
    dfig(13);clf;
    C = reshape(E,size(x));
    surf(p/pi, t/pi,double(C), 'EdgeColor', 'none');axis equal;view(0,-90);
    xlabel('phi');ylabel('theta');
    dfig(1);
    
    dfig(11);clf;
    C = reshape(E_test,size(x));
    surf(p/pi, t/pi,double(C), 'EdgeColor', 'none');axis equal;view(0,-90);
    xlabel('phi');ylabel('theta');
    dfig(1);
end
E = energy;
%%%%%%%%%%%%%%%%%%%%%%%%%%
