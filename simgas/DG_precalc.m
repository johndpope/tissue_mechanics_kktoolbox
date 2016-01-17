function [X_o, B, Y_ge] = DG_precalc(A_o, L_max, L_max_gef, gdim)

X_o = shp_sphere_gen(L_max,sqrt(A_o/4/pi));

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

[Y_ge]                   = ylk_cos_sin_bosh(p, t, L_max_gef);

wp = wp(:);wt = wt(:);w = wp.*wt;

B.Y = Y;
B.Y_P = Y_P;
B.Y_PP = Y_PP;
B.Y_T = Y_T;
B.Y_TT = Y_TT;
B.Y_TP = Y_TP;
B.p = p;
B.t = t;
B.w = w;