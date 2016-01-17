function [A, V,h,T, dA, E, H, wp, wt, p, t]  = plot_sh_notri(X_o, gdim)
% Calculates the shape properties of a shape defined using the spherical harmonics parameterization
% shape properties include:
%       Area, bending energy, volume, total curvature and Gaussian
%       curvature
% Author: Khaled Khairy 2006
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

%%
% % for ix = 1:size(Y_TT,3)
% % dfig(ix);
% % disp(ix);
% % clf;surf(p,t,Y(:,:,ix), 'EdgeColor','none');view(3);
% % end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdimp = length(p);gdimt = length(t);
c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
x = sum(c.*Y,3);
xp =(sum(c.*Y_P,3));
xt = (sum(c.*Y_T,3));
xpp = (sum(c.*Y_PP,3));xtt = (sum(c.*Y_TT,3));xtp = (sum(c.*Y_TP,3));

c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
y = sum(c.*Y,3);
yp =(sum(c.*Y_P,3));
yt = (sum(c.*Y_T,3));
ypp = (sum(c.*Y_PP,3));ytt = (sum(c.*Y_TT,3));ytp = (sum(c.*Y_TP,3));

c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
z = sum(c.*Y,3);
zp =(sum(c.*Y_P,3));
zt = (sum(c.*Y_T,3));
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


% % E = dot(Xt,Xt,2);
% % F = dot(Xt,Xp,2);
% % G = dot(Xp,Xp,2);
% % 
% % L = - dot(Xt,nt,2);
% % M = 1/2*(dot(Xt,np,2)+dot(Xp,nt,2));
% % N = - dot(Xp,np,2);
% % 
% % E = reshape(E, size(wp));
% % F = reshape(F, size(wp));
% % G = reshape(G, size(wp));
% % 
% % L = reshape(L, size(wp));
% % M = reshape(M, size(wp));
% % N = reshape(N, size(wp));
% % 
% % H = (E.*N + G.*L - 2*F.*M)./(2 * E.*G-F.*F);

%H = ((conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(-conj(xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))+(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp).*(-conj(xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp))))-(2.*conj(xt).*xp+2.*conj(yt).*yp+2.*conj(zt).*zp).*(1./2.*conj(xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))))./(2.*(conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp)-2.*(conj(xt).*xp+conj(yt).*yp+conj(zt).*zp).^2);

H = (( (xt).*xt+ (yt).*yt+ (zt).*zt).*(- (xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))+( (xp).*xp+ (yp).*yp+ (zp).*zp).*(- (xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp))))-(2.* (xt).*xp+2.* (yt).*yp+2.* (zt).*zp).*(1./2.* (xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))))./(2.*( (xt).*xt+ (yt).*yt+ (zt).*zt).*( (xp).*xp+ (yp).*yp+ (zp).*zp)-2.*( (xt).*xp+ (yt).*yp+ (zt).*zp).^2);

%% fix H spurious values
% % [tt wtt]                  = gaussquad(gdimt, 0, pi);
% % r = 10;
% % for ix = 1:r-1,
% % H(t==(tt(ix)))= H(t==tt(r));
% % H(t==(tt(end-ix + 1)))= H(t==tt(end-r));
% % end

%% %
h = 1./A*sum(sum(H.*wp.*wt.*SS));
K = ((-conj(xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))-conj(zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).*(-conj(xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))-conj(zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))-(1./2.*conj(xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.*conj(xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.*conj(zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))).^2)./((conj(xt).*xt+conj(yt).*yt+conj(zt).*zt).*(conj(xp).*xp+conj(yp).*yp+conj(zp).*zp)-(conj(xt).*xp+conj(yt).*yp+conj(zt).*zp).^2);
T = sum(sum(K.*wp.*wt.*SS))/4/pi;       % total curvature (constant for topology)
dA= sum(sum(wp.*wt.*2.*H.*SS));         % assuming unit distance D (since D goes away in some expressions below)
% Curvedness = (2.*H.^2-K).^(0.5);		% curvedness
fac = (H.^2-K).^(1/2);k_min = H-fac;k_max = H+fac;
E = 1/2*sum(sum((2.*H).^2.*wp.*wt.*SS))./8/pi;   % reduced bending energy
if plotflag,
    if sign(1/3.*sum(sum(wp.*wt.*v.*SS)))==-1, H = -H;end
%         minH = -0.8 ;  H(H<minH) = minH;
%         maxH = 0.2 ;  H(H>maxH) = maxH;
    surf(x,y,z,double(H), 'EdgeColor', 'none');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
