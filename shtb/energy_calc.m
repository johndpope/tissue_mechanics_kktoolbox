function E = energy_calc(X_o)
% Calculates the energy by quadrature
% Author: Khaled Khairy 2010
%%% XXX Works as integration but very slow, the results are still inaccurate --- use as check only XXX
[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);
tol = 1e-2;
E = dblquad(@energy_integrnd_inner,0,2*pi, 0, pi);

%     function dE1 = energy_integrnd_outer(p)
%         
%         dE1 = quadgk(@energy_integrnd_inner, 0, pi);
        function dE = energy_integrnd_inner(p,t)
            p = p(:)';
            t = t(:)';
            if t>-1e-6 & t<1e-6, t = 1e-6;end
            if t<pi+1e-6 & t>pi-1e-6, t = pi-1e-6;end
            
            [Y P]                   = ylk_cos_sin_bosh(p, t, L_max);
            Y_P                     = ylk_cos_sin_dphi_bosh(p, t, L_max, P);
            [Y_T P_T]               = ylk_cos_sin_dtheta_bosh(p, t, L_max, P);
            Y_PP 					= ylk_cos_sin_dphiphi_bosh(p, t, L_max, P);
            Y_TT 					= ylk_cos_sin_dthetatheta_bosh(p, t, L_max,P_T);
            Y_TP 					= ylk_cos_sin_dthetaphi_bosh(p, t, L_max, P_T);
            
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
            
            %         X  =[x(:) y(:) z(:)];
            Xt = [xt(:) yt(:) zt(:)];
            Xp = [xp(:) yp(:) zp(:)];
            %         Xpp = [xpp(:) ypp(:) zpp(:)];
            %         Xtp = [xtp(:) ytp(:) ztp(:)];
            %         Xtt = [xtt(:) ytt(:) ztt(:)];
            
            SS1 = (cross(Xt,Xp,2));
            SS = sqrt(SS1(:,1).^2 + SS1(:,2).^2+SS1(:,3).^2);
            
            H = (( (xt).*xt+ (yt).*yt+ (zt).*zt).*(- (xp).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (yp).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))- (zp).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp))))+( (xp).*xp+ (yp).*yp+ (zp).*zp).*(- (xt).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (yt).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))- (zt).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp))))-(2.* (xt).*xp+2.* (yt).*yp+2.* (zt).*zp).*(1./2.* (xt).*((ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (yt).*((ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (zt).*((xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytp.*zp+yt.*zpp-ztp.*yp-zt.*ypp)+2.*(zt.*xp-xt.*zp).*(ztp.*xp+zt.*xpp-xtp.*zp-xt.*zpp)+2.*(xt.*yp-yt.*xp).*(xtp.*yp+xt.*ypp-ytp.*xp-yt.*xpp)))+1./2.* (xp).*((ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(yt.*zp-zt.*yp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (yp).*((ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(zt.*xp-xt.*zp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))+1./2.* (zp).*((xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(1./2)-1./2.*(xt.*yp-yt.*xp)./((yt.*zp-zt.*yp).^2+(zt.*xp-xt.*zp).^2+(xt.*yp-yt.*xp).^2).^(3./2).*(2.*(yt.*zp-zt.*yp).*(ytt.*zp+yt.*ztp-ztt.*yp-zt.*ytp)+2.*(zt.*xp-xt.*zp).*(ztt.*xp+zt.*xtp-xtt.*zp-xt.*ztp)+2.*(xt.*yp-yt.*xp).*(xtt.*yp+xt.*ytp-ytt.*xp-yt.*xtp)))))./(2.*( (xt).*xt+ (yt).*yt+ (zt).*zt).*( (xp).*xp+ (yp).*yp+ (zp).*zp)-2.*( (xt).*xp+ (yt).*yp+ (zt).*zp).^2);
            dE = 1/2*(2.*H).^2.*SS./8/pi;   % reduced bending energy
            if isnan(dE),warning('check');end
        end
%    end
end