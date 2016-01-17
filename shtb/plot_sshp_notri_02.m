function [E, A, V,h,T, dA]  = plot_sshp_notri(X_o, gdim)
% Calculates the shape properties of a shape defined using the spherical harmonics parameterization
% and plots it.
% Note: this is compatible with the symbolic construction of the basis,
% i.e. the "strict" version
% Author: Khaled A. Khairy 2010
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1;plotflag = 1;
if nargin == 1,
    gdim = 60;
end
%% generate the symbolic basis
syms PI t p positive
syms x
L_max = get_L_max(X_o);
[L K] = indices_gen(1:(L_max + 1)^2);
disp('Generating basis expressions...')
for ix = 1:length(L),           % loop over all pairs of legal L and K values
    str = sprintf('syms xc%d yc%d zc%d real',ix,ix,ix);eval(str);
    ylk(ix) = sym_YLK(L(ix),K(ix));
    ylk_p(ix) = diff(ylk(ix),p);
    ylk_pp(ix) = diff(ylk_p(ix),p);
    ylk_t(ix) = diff(ylk(ix), t);
    ylk_tp(ix) = diff(ylk_t(ix), p);
    ylk_tt(ix) = diff(ylk_t(ix), t);
end
disp('Done');
%% generate the general surface parameterization in functional form S
Sx = 0;Sy = 0;Sz = 0;
for ix = 1:length(L),
    str = sprintf('Sx = Sx + xc%d * ylk(ix);',ix);eval(str);
    str = sprintf('Sy = Sy + yc%d * ylk(ix);',ix);eval(str);
    str = sprintf('Sz = Sz + zc%d * ylk(ix);',ix);eval(str);
end
S = [Sx Sy Sz];
%% obtain the derivatives of the surface function w.r.t. theta and phi
St = diff(S,t);
Sp = diff(S,p);

%% the surface differentials, surface normal and area element still general
SS1 = cross(St,Sp,2);
SS = sqrt(SS1(:,1).^2 + SS1(:,2).^2+SS1(:,3).^2);
n = SS1./SS(:,ones(1,3));
v  = dot(S,n,2);    % volume element
%% Perform Gaussian quadrature to obtain final values
%% generate the Gaussian quadrature angles and weights
[t wt]                  = gaussquad(gdim, 0, pi);
[p wp]                  = gaussquad(gdim,0,2*pi);
[p t]                   = meshgrid(p,t);
[wp wt]                 = meshgrid(wp, wt);
t = t(:);p = p(:);
%% assign the coefficients their values
[xclks yclks zclks] = get_xyz_clks(X_o);
for ix = 1:length(L),
    str = sprintf('xc%d = %.6f;yc%d = %.6f;zc%d = %.6f;',ix,xclks(ix),ix,yclks(ix), ix, zclks(ix));eval(str);
end
%% evaluate the basis expressions and the derivatives
PI = pi;
SS1 = eval(vectorize(SS1));
SS  = eval(vectorize(SS));
n   = eval(vectorize(SS));

%%  Calculate the properties     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if sign(1/3.*sum(sum(wp.*wt.*v.*SS)))==-1, H = -H;end
    minH = -3.5736 ;  H(H<minH) = minH;
    maxH = 2.3223 ;  H(H>maxH) = maxH;
        surf(x,y,z,double(H), 'EdgeColor', 'none');axis equal;axis tight; lighting gouraud;light;end












