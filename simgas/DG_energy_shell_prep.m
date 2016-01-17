function [phys] = DG_energy_shell_prep(X_o,phys)
%% process input
Y = phys.B.Y;
Y_P = phys.B.Y_P;
Y_T = phys.B.Y_T;
Y_PP = phys.B.Y_PP;
Y_TT = phys.B.Y_TT;
Y_TP = phys.B.Y_TP;
p    = phys.B.p;
t    = phys.B.t;
w    = phys.B.w;
X_o = recover_parms(X_o,phys);
X_o = over_Basis(X_o);
[xclks yclks zclks] = get_xyz_clks(X_o);
%% generate x y z, and their individual derivatives with respect to theta and phi
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
%% generate X surface vector and derivatives of X
X  =[x(:) y(:) z(:)];
Xt = [xt(:) yt(:) zt(:)];       % m_l1
Xp = [xp(:) yp(:) zp(:)];       % m_l2
Xpp = [xpp(:) ypp(:) zpp(:)];
Xtp = [xtp(:) ytp(:) ztp(:)];
Xtt = [xtt(:) ytt(:) ztt(:)];

%% generate unit tangent vectors
Xtnorm = sqrt(kk_dot(Xt,Xt));
Xpnorm = sqrt(kk_dot(Xp,Xp));

tt = Xt./Xtnorm(:,ones(1,3));
tp = Xp./Xpnorm(:,ones(1,3));
%% calculate coefficients of the 1st fundamental form (components of g_lij)
E = kk_dot(Xt,Xt);          % g_l11
F = kk_dot(Xt,Xp);          % g_l12 and g_l21
G = kk_dot(Xp,Xp);          % g_l22
%% calculate differential area element and surface normal
SS = (kk_cross(Xt,Xp));
SSn = sqrt(E.*G-F.*F);
n = SS./SSn(:,ones(1,3));       % m_l3
%% set up covariant basis (with letter l), and contravariant basis with the letter u
m_l1 = Xt;
m_l2 = Xp;
m_l3 = n;
bb = kk_dot(m_l1, kk_cross(m_l2,m_l3));bb = bb(:);  % formulation as in AMoS p.658
bb = kk_norm(kk_cross(m_l1,m_l2));                  % bb is essentially SSn
bb = bb(:,ones(1,3));
m_u1 = kk_cross(Xp,n)./bb;
m_u2 = kk_cross(n,Xt)./bb;
m_u3 = n;
%% set up the metric tensor components
% covariant components
g_l11 = E;
g_l22 = G;
g_l12 = F;
g_l13 = kk_dot(m_l1, m_l3);
g_l23 = kk_dot(m_l2, m_l3);
g_l33 = kk_dot(m_l3, m_l3);
% contravariant metric tensor components g_uij
g_u11 = kk_dot(m_u1, m_u1);
g_u12 = kk_dot(m_u1, m_u2);
g_u22 = kk_dot(m_u2, m_u2);
g_u33 = kk_dot(m_u3, m_u3);
g_u13 = kk_dot(m_u1, m_u3);
g_u23 = kk_dot(m_u2, m_u3);
%% calculate the coefficients of the second fundamental form
L = kk_dot(Xtt,n);
M = kk_dot(Xtp,n);
N = kk_dot(Xpp,n);
SO = [(G.*L-F.*M)./(E.*G-F.*F) (E.*M-F.*L)./(E.*G-F.*F) (G.*M-F.*N)./(E.*G-F.*F) (E.*N-F.*M)./(E.*G-F.*F)]; % the shape operator
%SO = [(G.*L-F.*M) (E.*M-F.*L) (G.*M-F.*N) (E.*N-F.*M)]; % the shape operator

% fac = 1;%SSn.^2;
k_u1_lo1 = SO(:,1);
k_u1_lo2 = SO(:,3);
k_u2_lo1 = SO(:,2);
k_u2_lo2 = SO(:,4);

% k_lo11 = kk_dot(m_l3,Xtt);
% k_lo12 = kk_dot(m_l3,Xtp);
% k_lo21 = kk_dot(m_l3,Xtp);
% k_lo22 = kk_dot(m_l3,Xpp);
% 
%  
% k_u1_lo1 = g_u11.*k_lo11 + g_u12.*k_lo21;
% k_u1_lo2 = g_u11.*k_lo12 + g_u12.*k_lo22;
% k_u2_lo1 = g_u12.*k_lo11 + g_u22.*k_lo21;
% k_u2_lo2 = g_u12.*k_lo12 + g_u22.*k_lo22;


%% calculate  geometric properties: Area, Volume, reduced volume, Curvatures and self-check total Gaussian curvature
V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
A = sum(w.*SSn);
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
K = (L.*N - M.*M)./(E.*G-F.*F);
h = 1./A.*sum(H(:).*w.*SSn);
T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
%% store preliminary quantities in phys
    phys.A = A;
    phys.V = V;
    phys.C_o = H;
    phys.Xt_o = Xt;
    phys.Xp_o = Xp;
    phys.E  = E;
    phys.G  = G;
    phys.F  = F;
    phys.L = L;
    phys.M = M;
    phys.N = N;
    phys.SS_o = SS;
    phys.SSn_o = SSn;
    phys.m_l1 = m_l1;
    phys.m_l2 = m_l2;
    phys.m_l3 = m_l3;
    phys.m_u1 = m_u1;
    phys.m_u2 = m_u2;
    phys.m_u3 = m_u3;
    
    phys.g_u11 = g_u11;
    phys.g_u12 = g_u12;
    phys.g_u13 = g_u13;
    phys.g_u22 = g_u22;
    phys.g_u23 = g_u23;
    phys.g_u33 = g_u33;
    
    phys.g_l11 = E;
    phys.g_l12 = F;
    phys.g_l13 = g_l13;
    phys.g_l22 = G;
    phys.g_l23 = g_l23;
    phys.g_l33 = g_l33;
    
    phys.tt = tt;
    phys.tp = tp;
    
    phys.k_lo11 = SO(:,1);
    phys.k_lo12 = SO(:,3);
    phys.k_lo21 = SO(:,2);
    phys.k_lo22 = SO(:,4);

    phys.k_u1_lo1 = k_u1_lo1;
    phys.k_u1_lo2 = k_u1_lo2;
    phys.k_u2_lo1 = k_u2_lo1;
    phys.k_u2_lo2 = k_u2_lo2;
%%%%%%%%%%%%%%%%%%%%%    
function r = kk_norm(a)
%%%% a is a nx3 vector
%%%% r is a nx1 vector of the norm of every row in a
r = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);

