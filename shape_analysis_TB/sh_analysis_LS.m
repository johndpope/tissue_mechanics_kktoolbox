function [clks] = sh_analysis_LS(X, L_max, flag)

%%%%%%%%%%%%%%%%%%%%%% Convert to spherical polar coordinates
x = X(:,1)'; y = X(:,2)'; z = X(:,3)';
R = sqrt(x.^2+y.^2+z.^2);
phi = atan2(y,x);
theta = atan2(sqrt(x.^2+y.^2),z);

%%%%%%%%%%%%%%%%%%%%%% Solve Linear Least squares
warning off MATLAB:divideByZero;
[L, K ] = indices_gen(1:(L_max + 1)^2);
%2 * length(L);
M = length(L); %% number of functions in expansion
N = length(R); %% number of data points
A  = zeros(N, M);
for S = 1:length(L)
    A(:, S) = ylk_cos_sin_old(L(S),K(S),phi, theta)';
end
clks = inv(A'*A)*(A'*R');

% [U, S, V] = svd(A, 0);
% invS = 1./(S);
% invS(invS==inf) = 0;
% clks = (V*invS) * (U'*R');

%%%%%%%%%%%%%%%%%%%%% Calculate initial quantities
if flag
    gdim = 120;
    global Y_LK P_LK Y_LK_phi Y_LK_theta P_LK_T Y_PP Y_TT Y_TP wp wt
    [t wt]                  = gaussquad(gdim, 0, pi);
    [p wp]                  = gaussquad(gdim,0,2*pi);
    [p t]                   = meshgrid(p,t);
    [wp wt]                 = meshgrid(wp, wt);
    [Y_LK P_LK]				= precalc_ylk_cos_sin(p, t, L_max);
    Y_LK_phi 				= precalc_ylk_cos_sin_dphi(p, t, L_max);
    [Y_LK_theta P_LK_T] 	= precalc_ylk_cos_sin_dtheta(p, t, L_max);
    Y_PP 					= precalc_ylk_cos_sin_dphiphi(p, t, L_max);
    Y_TT 					= precalc_ylk_cos_sin_dthetatheta(p, t, L_max);
    Y_TP 					= precalc_ylk_cos_sin_dthetaphi(p, t, L_max);

    st  = sin(t);sp = sin(p); ct = cos(t); cp = cos(p);
    ctsq = ct.^2; cpsq = cp.^2; stsq = st.^2; spsq = sp.^2;
    ctcu = ct.^3; ctqad = ct.^4;
    %%%%%%%%%%%%%%%%%%%%%%% Calculate shape properties
    x = clks';h = length(x);x = x(ones(gdim^2,1),:);x = reshape(x,gdim, gdim, h); % expensive
    r_o = 0;
    R = r_o + sum(x.*Y_LK,3);Rp =(sum(x.*Y_LK_phi,3));Rt = (sum(x.*Y_LK_theta,3));
    Rpp = (sum(x.*Y_PP,3));Rtt = (sum(x.*Y_TT,3));Rtp = (sum(x.*Y_TP,3));
    deltaR = Rtt + ct./st.*Rt + Rpp./stsq;
    divRsq = Rt.^2 + Rp.^2./stsq;
    term3 = Rt.^2.*Rtt + Rt.*Rp.*Rtp./stsq -Rt.*Rp.^2 .* ct./st.^3 + Rt.*Rp.*Rtp./stsq + Rp.^2 .* Rpp./stsq.^2;
    integrand1 = 2.*R -deltaR + (R.*divRsq + term3)./(R.^2 + divRsq);
    H = integrand1/2;
    H_min = -20; H_max = 20;
    H(H>H_max) = H_max; H(H<H_min) = H_min;
    %%%%%%%%%%%%%%%%%%%%%%% Look at the shape corresponding to c
    r = R.*st;x = r.*cp;y = r.*sp;z = R.*ct;
    surf(x,y,z,(Rp.^2 + Rt.^2));grid off;daspect([1 1 1]);hold off;colorbar;rotate3d;drawnow;

end