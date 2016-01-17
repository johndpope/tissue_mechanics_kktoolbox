%% % test symbolically the second term of DGT on p 665 AMoS
clc;clear all;
syms H E G F L N M real
syms H_ E_ G_ F_ L_ M_ real
syms N__ g g_ real
% g = sqrt(E*G-F^2);
k = [G*L-F*M G*M-F*N;E*M-F*L E*N-F*M]./g^2;
H = trace(k);

% g_ = sqrt(E_*G_-F_^2);
k_ = [G_*L_-F_*M_ G_*M_-F_*N__;E_*M_-F_*L_ E_*N__-F_*M_]./g_^2;
H_ = trace(k_);% H_ = E_*N_+G_*L_-2*F_*M_/g_^2;

%% define tangent vectors 
syms m1 m2 m3 m_1 m_2 m_3 mu1 mu2 mu3 m_u1 m_u2 m_u3 xt yt zt xp yp zp t p real
syms xt_ yt_ zt_ xp_ yp_ zp_ real
syms cm1m2 cm2m3 cm3m1 cm_1m_2 cm_2m_3 cm_3m_1 
m1 = [xt yt zt];
m2 = [xp yp zp];
m3 = cross(m1,m2)/g;
mu1 = cross(m2,m3)/g;
mu2 = cross(m3,m1)/g;
mu3 = m3;

m_1 = [xt_ yt_ zt_];
m_2 = [xp_ yp_ zp_];
m_3 = cross(m_1, m_2)/g_;
m_u1 = cross(m_2,m_3)/g_;
m_u2 = cross(m_3,m_1)/g_;
m_u3 = m_3;

E = dot(m1,m1);
G = dot(m2,m2);
F = dot(m1,m2);

E_ = dot(m_1,m_1);
G_ = dot(m_2,m_2);
F_ = dot(m_1,m_2);


%% metric tensor
g_11 = E_;
g_22 = G_;
g_12 = F_;
% contravariant metric tensor components g_uij
g_u11 = dot(m_u1, m_u1);
g_u12 = dot(m_u1, m_u2);
g_u22 = dot(m_u2, m_u2);

g11 = E;
g22 = G;
g12 = F;
% contravariant metric tensor components guij
gu11 = dot(mu1, mu1);
gu12 = dot(mu1, mu2);
gu22 = dot(mu2, mu2);
%% Calculate the deformed and undeformed mixed components of the shape tensor kappa
k_lo11 = k_(1,1);
k_lo12 = k_(1,2);
k_lo21 = k_(2,1);
k_lo22 = k_(2,2);
k_u1_lo1 = g_u11.*k_lo11 + g_u12.*k_lo21;
k_u1_lo2 = g_u11.*k_lo12 + g_u12.*k_lo22;
k_u2_lo1 = g_u12.*k_lo11 + g_u22.*k_lo21;
k_u2_lo2 = g_u12.*k_lo12 + g_u22.*k_lo22;

klo11 = k(1,1);
klo12 = k(1,2);
klo21 = k(2,1);
klo22 = k(2,2);
ku1_lo1 = gu11.*klo11 + gu12.*klo21;
ku1_lo2 = gu11.*klo12 + gu12.*klo22;
ku2_lo1 = gu12.*klo11 + gu22.*klo21;
ku2_lo2 = gu12.*klo12 + gu22.*klo22;

d11 = (ku1_lo1 - k_u1_lo1);
d12 = (ku1_lo2 - k_u1_lo2);
d21 = (ku2_lo1 - k_u2_lo1);
d22 = (ku2_lo2 - k_u2_lo2);
DGT2 = (...
    d11.*kron(m1',m_u1) + ...
    d12.*kron(m1',m_u2) + ...
    d21.*kron(m2',m_u1) + ...
    d22.*kron(m2',m_u2) );
BC = DGT2'*DGT2;%kk_mx_mult(kk_transpose(DGT2), DGT2);

res = sum(BC(:))/(H-H_).^2;





