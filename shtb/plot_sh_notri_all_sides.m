function [E, A, V,h,T, dA]  = plot_sh_notri_all_sides(X_o, gdim)
% Calculates the shape properties of a shape defined using the spherical harmonics parameterization
% shape properties include: 
%       Area, bending energy, volume, total curvature and Gaussian
%       curvature
% Will initialize itself if precalculated (global) quantities are absent.
% Author: Khaled Khairy 2006
warning off; global Y_LK P_LK Y_LK_phi Y_LK_theta P_LK_T Y_PP Y_TT Y_TP gdimp gdimt wp wt p t;warning on;
verbose = 1;plotflag = 1;
if nargin == 1,
    gdim = 60;
end
nc = length(X_o)/3;
xclks = X_o(1:nc);
yclks = X_o(nc+1:2*nc);
zclks = X_o(2*nc+1:end);

verbose = 1;
L_max = round(sqrt(length(xclks))-1);
if size(Y_LK_theta,1)~=gdim || size(Y_LK_theta,3)~=(L_max+1)^2,
    [t wt]                  = gaussquad(gdim, 0, pi);
    [p wp]                  = gaussquad(gdim,0,2*pi);
    [p t]                   = meshgrid(p,t);
    [wp wt]                 = meshgrid(wp, wt);
    if verbose, disp('Generating Y_LKs');end
    [Y_LK P_LK]				= precalc_ylk_cos_sin(p, t, L_max); 
    if verbose, disp('Generating Y_LK_phi');end
    Y_LK_phi 				= precalc_ylk_cos_sin_dphi(p, t, L_max);
    if verbose, disp('Generating Y_LK_theta');end
    [Y_LK_theta P_LK_T] 	= precalc_ylk_cos_sin_dtheta(p, t, L_max);
    if verbose, disp('Generating Y_LK_pp');end
    Y_PP 					= precalc_ylk_cos_sin_dphiphi(p, t, L_max);
    if verbose, disp('Generating Y_LK_tt');end
    Y_TT 					= precalc_ylk_cos_sin_dthetatheta(p, t, L_max);
    if verbose, disp('Generating Y_LK_tp');end
    Y_TP 					= precalc_ylk_cos_sin_dthetaphi(p, t, L_max);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdimp = length(p);gdimt = length(t);
c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
x = sum(c.*Y_LK,3);
xp =(sum(c.*Y_LK_phi,3));
xt = (sum(c.*Y_LK_theta,3));
xpp = (sum(c.*Y_PP,3));xtt = (sum(c.*Y_TT,3));xtp = (sum(c.*Y_TP,3));

c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
y = sum(c.*Y_LK,3);
yp =(sum(c.*Y_LK_phi,3));
yt = (sum(c.*Y_LK_theta,3));
ypp = (sum(c.*Y_PP,3));ytt = (sum(c.*Y_TT,3));ytp = (sum(c.*Y_TP,3));

c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
z = sum(c.*Y_LK,3);
zp =(sum(c.*Y_LK_phi,3));
zt = (sum(c.*Y_LK_theta,3));
zpp = (sum(c.*Y_PP,3));ztt = (sum(c.*Y_TT,3));ztp = (sum(c.*Y_TP,3));

%%%%%%%%%%%%%%%%%%%%%%%% Calculate the geometrically important quantities

X  =[x(:) y(:) z(:)];
Xt = [xt(:) yt(:) zt(:)];
Xp = [xp(:) yp(:) zp(:)];
Xpp = [xpp(:) ypp(:) zpp(:)];
Xtp = [xtp(:) ytp(:) ztp(:)];
Xtt = [xtt(:) ytt(:) ztt(:)];

SS1 = (cross(Xt,Xp,2));
SS = sqrt(SS1(:,1).^2 + SS1(:,2).^2+SS1(:,3).^2);
n = SS1./SS(:,ones(1,3));
v  = dot(X,n,2);
v = reshape(v,size(wp));
SS = reshape(sqrt(SS1(:,1).^2 + SS1(:,2).^2+SS1(:,3).^2),size(wp));
V = abs(1/3.*sum(sum(wp.*wt.*v.*SS)));
A = sum(sum(wp.*wt.*SS));

H = ((conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(-conj(xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))+(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp).*(-conj(xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp))))-(2.*conj(xt).*xp+2.*conj(yt).*yp+2.*conj(zt).*zp).*(1./2.*conj(xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))))./(2.*(conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp)-2.*(conj(xt).*xp+conj(yt).*yp+conj(zt).*zp).^2);
h = 1./A*sum(sum(H.*wp.*wt.*SS));
K = ((-conj(xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).*(-conj(xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))-(1./2.*conj(xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).^2)./((conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp)-(conj(xt).*xp+conj(yt).*yp+conj(zt).*zp).^2);
T = sum(sum(K.*wp.*wt.*SS))/4/pi;       % total curvature (constant for topology)
dA= sum(sum(wp.*wt.*2.*H.*SS));         % assuming unit distance D (since D goes away in some expressions below)
% Curvedness = (2.*H.^2-K).^(0.5);		% curvedness
% fac = (H.^2-K).^(1/2);k_min = H-fac;k_max = H+fac;

%%%%%%%%%%%%%%%%% Calculate the Energies
E = 1/2*sum(sum((2.*H).^2.*wp.*wt.*SS))./8/pi;   % reduced bending energy


%%%%%%%%%%%%%%%% Display results
if verbose,
    v_red = V./(4/3*pi*(A/4/pi).^(3/2));
    str = sprintf('Area: %.3f \nVolume: %.3f \nRedVol: %.3f \nBending E: %.3f \nh: %.3f \nT:%.3f',...
        A, V,v_red,E, h,T);disp(str);end

if plotflag, 
    clf;
    subplot(2,2,1);surf(x,y,z,double(H), 'EdgeColor', 'none');axis equal;axis tight; lighting gouraud;light;view(0,0);xlabel('x');ylabel('y');
    subplot(2,2,2);surf(x,y,z,double(H), 'EdgeColor', 'none');axis equal;axis tight; lighting gouraud;light;view(90,0);xlabel('x');ylabel('y');
    subplot(2,2,3);surf(x,y,z,double(H), 'EdgeColor', 'none');axis equal;axis tight; lighting gouraud;light;view(0,90);xlabel('x');ylabel('y');
    subplot(2,2,4);surf(x,y,z,double(H), 'EdgeColor', 'none');axis equal;axis tight; lighting gouraud;light;view(3);xlabel('x');ylabel('y');
end













