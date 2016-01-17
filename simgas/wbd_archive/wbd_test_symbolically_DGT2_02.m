%% % test symbolically the second term of DGT on p 665 AMoS
clc;clear all;
%% define the tangent vectors
syms xt yt zt xtt ytt ztt real
syms xp yp zp xpp ypp zpp real
syms xtp ytp ztp g real
Xtt = [xtt ytt ztt];
Xpp = [xpp ypp zpp];
Xtp = [xtp ytp ztp];
m1 = [xt yt zt];
m2 = [xp yp zp];
E = dot(m1,m1);
G = dot(m2,m2);
F = dot(m1,m2);

% g = sqrt(E*G-F^2);

m3 = cross(m1,m2)/g;
n = m3;

L = dot(Xtt,n);
M = dot(Xtp,n);
N = dot(Xpp,n);

mu1 = cross(m2, m3)/g;
mu2 = cross(m3, m1)/g;
%% define shape operator and the local mean curvature
k = [G*L-F*M G*M-F*N;E*M-F*L E*N-F*M]./g^2;
H = trace(k);
Eb1 = 2*(H)^2;
%% metric tensor
g11 = E;
g22 = G;
g12 = F;
% contravariant metric tensor components guij
gu11 = dot(mu1, mu1);
gu12 = dot(mu1, mu2);
gu22 = dot(mu2, mu2);
%% Calculate the mixed components of the shape tensor kappa
klo11 = k(1,1);
klo12 = k(1,2);
klo21 = k(2,1);
klo22 = k(2,2);
ku1_lo1 = gu11.*klo11 + gu12.*klo21;
ku1_lo2 = gu11.*klo12 + gu12.*klo22;
ku2_lo1 = gu12.*klo11 + gu22.*klo21;
ku2_lo2 = gu12.*klo12 + gu22.*klo22;

d11 = (ku1_lo1 - 0);
d12 = (ku1_lo2 - 0);
d21 = (ku2_lo1 - 0);
d22 = (ku2_lo2 - 0);
% DGT2 = (...
%     d11.*kron(m1',m_u1) + ...
%     d12.*kron(m1',m_u2) + ...
%     d21.*kron(m2',m_u1) + ...
%     d22.*kron(m2',m_u2) );
DGT2 = [d11 0;0 d22];
BC = DGT2'*DGT2;%kk_mx_mult(kk_transpose(DGT2), DGT2);

Eb2 = sum(BC(:));
res = Eb1/Eb2;
disp(res)




