function [E, A, V, v_red, AreaDiff, Eb, EADE, E_stretch, E_shear, T, H, phys] = ...
    DG_energy_shell_08_DGT(X_o,phys)
%% process input
Eb = 0;
EADE = 0;
E_stretch = 0;
E_shear = 0;
AreaDiff = 0;

C_o = phys.C_o;

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
bb = kk_norm(kk_cross(m_l1,m_l2));                  % bb is essentially SSn
bb = bb(:,ones(1,3));
m_l3 = kk_cross(Xt,Xp)./bb;

% bb = kk_dot(m_l1, kk_cross(m_l2,m_l3));bb = bb(:);  % formulation as in AMoS p.658
m_u1 = kk_cross(Xp,m_l3)./bb;
m_u2 = kk_cross(m_l3,Xt)./bb;
m_u3 = n;
%% generate unit tangent vectors -- for treatment as in AoPaS
Xtnorm = sqrt(kk_dot(Xt,Xt));
Xpnorm = sqrt(kk_dot(Xp,Xp));

tt = Xt./Xtnorm(:,ones(1,3));
tp = Xp./Xpnorm(:,ones(1,3));
tn = kk_cross(Xt,Xp)./Xpnorm(:,ones(1,3))./Xtnorm(:,ones(1,3));
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
% g_u33 = kk_dot(m_u3, m_u3);
% g_u13 = kk_dot(m_u1, m_u3);
% g_u23 = kk_dot(m_u2, m_u3);
%% calculate the coefficients of the second fundamental form
L = kk_dot(Xtt,n);
M = kk_dot(Xtp,n);
N = kk_dot(Xpp,n);
SO = [(G.*L-F.*M)./(E.*G-F.*F) (E.*M-F.*L)./(E.*G-F.*F) (G.*M-F.*N)./(E.*G-F.*F) (E.*N-F.*M)./(E.*G-F.*F)]; % the shape operator
% SO = [(G.*L-F.*M) (E.*M-F.*L) (G.*M-F.*N) (E.*N-F.*M)]; % the shape operator
%% calculate  geometric properties: Area, Volume, reduced volume, Curvatures and self-check total Gaussian curvature
V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
A = sum(w.*SSn);
Vo = 4/3*pi*(A/4/pi)^(3/2);
v_red = V/Vo;
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
K = (L.*N - M.*M)./(E.*G-F.*F);
h = 1./A.*sum(H(:).*w.*SSn);
T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
%% DGT 1: Cauchy-Green deformation tensor -> principal stretches
%%% first term corresponds to in-plane deformations (stretch and shear)
%%% vectorized version
DGT1 = kk_kron(m_l1, phys.m_u1) + ...
      kk_kron(m_l2, phys.m_u2);
%   FT  = kk_transpose(DGT1);
% C = kk_mx_mult(FT,DGT1);% right Cauchy-Green deformation tensor
% [res] = real(eig3(reshape(C',3,3,size(C,1)))); % uses vectorized Cardan's formula for root finding
% res = sort(res',2);
% tol = 1e-4;indx = find(res(1,:)<(tol) & res(1,:)>(-tol));
% if isempty(indx), warning('empty indx');res = sort(res,2,'descend');res(:,3) = [];
% else res(:,indx(1)) = [];
% end
% lambda = sqrt(res); % take square root for finding the principal stretches
% beta  = real((lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2));
% ai = real(lambda(:,1).*lambda(:,2) - 1);     % area_invariant (as function of principal stretches
% 
% 
% E_stretch = phys.D * phys.k_stretch/2 * ...
%     sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;
% 
% E_shear = phys.D * phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% %% works: calculate bending energy as in RBC paper
% Eb = phys.D^3/2*sum((2.*H(:)-2.*C_o(:)).^2.*w.*SSn)./8/pi;   % reduced bending energy

%% TESTING DGT 2: Calculate term2 of F in AMoS p. 665
% fac = 1;%SSn.^2;
k_u1_lo1 = SO(:,1);
k_u1_lo2 = SO(:,3);
k_u2_lo1 = SO(:,2);
k_u2_lo2 = SO(:,4);
zer = zeros(length(Xt),1);
d11 = (k_u1_lo1 - phys.k_u1_lo1);
d12 = (k_u1_lo2 - phys.k_u1_lo2);
d21 = (k_u2_lo1 - phys.k_u2_lo1);
d22 = (k_u2_lo2 - phys.k_u2_lo2);
DGT2 = phys.D*[d11 zer zer zer d22 zer zer zer zer];% + ...
%     d12(:,ones(1,9)).*kk_kron(phys.m_l1,phys.m_u2) + ...
%     d21(:,ones(1,9)).*kk_kron(phys.m_l2,phys.m_u1);

%%% the equivalence to the bending energy term from rbc work
% BC = kk_mx_mult(kk_transpose(DGT2), DGT2);
% Eb2 =  sum(sum(BC,2).*w.*SSn)/8/pi;   % reduced bending energy
%% TESTING: DGT 3: third term of DGT
DGT3 = kk_kron(m_l3, phys.m_u3);%%% Calculate  DGT = kron(defgi, undefgi) 
%% TESTING: integration of the total energy
DGT = DGT1 + DGT2 + DGT3;
detF = kk_mx_det(DGT);
trFTF = kk_trace(kk_mx_mult(kk_transpose(DGT), DGT));
PHI = phys.lambda/2*log(detF.^2) + phys.miu/2.*(trFTF-3)-phys.miu/2.*log(detF); % Strain energy as in Conte et al.2007
PHI(isnan(PHI)) = 0;
PHI(isinf(PHI)) = 0;
E = phys.D*sum(PHI.*w.*SSn);     % perform the integration
% % E = Eb + E_stretch + E_shear;
%% --------------------------------------------------------------------
function r = kk_trace(a)
%%% a is nx9
%%% returns nx1
r = a(:,1)+a(:,5)+a(:,9);
function r = kk_transpose(a)
%%% a is nx9
%%% returns nx9
r = [a(:,1) a(:,4) a(:,7) a(:,2) a(:,5) a(:,8) a(:,3) a(:,6) a(:,9)];
function r = kk_mx_det(a)
%%% a is a nx9 array (each row represents components of a 2nd rank tensor
%%% returns nx1
r = a(:,1).*a(:,5).*a(:,9) +...
    a(:,4).*a(:,8).*a(:,3) +...
    a(:,7).*a(:,2).*a(:,6) -...
    a(:,7).*a(:,5).*a(:,3) -...
    a(:,4).*a(:,2).*a(:,9) -...
    a(:,1).*a(:,8).*a(:,6);
function r = kk_mx_mult(a,b)
%%% a and b are assumed to be nx9 arrays -- we are vectorizing matrix
%%% multiplication.
%%% returns an nx9 array
r1 = a(:,1).*b(:,1) + a(:,4).*b(:,2) + a(:,7).*b(:,3);
r4 = a(:,1).*b(:,4) + a(:,4).*b(:,5) + a(:,7).*b(:,6);
r7 = a(:,1).*b(:,7) + a(:,4).*b(:,8) + a(:,7).*b(:,9);
r2 = a(:,2).*b(:,1) + a(:,5).*b(:,2) + a(:,8).*b(:,3);
r5 = a(:,2).*b(:,4) + a(:,5).*b(:,5) + a(:,8).*b(:,6);
r8 = a(:,2).*b(:,7) + a(:,5).*b(:,8) + a(:,8).*b(:,9);
r3 = a(:,3).*b(:,1) + a(:,6).*b(:,2) + a(:,9).*b(:,3);
r6 = a(:,3).*b(:,4) + a(:,6).*b(:,5) + a(:,9).*b(:,6);
r9 = a(:,3).*b(:,7) + a(:,6).*b(:,8) + a(:,9).*b(:,9);
r = [r1 r2 r3 r4 r5 r6 r7 r8 r9];
function [r] = kk_kron(a,b)
%%% a and b are assumed to be nx3 arrays
%%% returns an nx9 array
r  = [a(:,1).*b(:,1) ...
    a(:,2).*b(:,1) ...
    a(:,3).*b(:,1) ...
    a(:,1).*b(:,2) ...
    a(:,2).*b(:,2) ...
    a(:,3).*b(:,2) ...
    a(:,1).*b(:,3) ...
    a(:,2).*b(:,3) ...
    a(:,3).*b(:,3)];
function [r na] = kk_norm(a)
%%%% a is a nx3 vector
%%%% r is a nx1 vector of the norm of every row in a
%%%% na = normalized a
r = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);
if nargout ==2, na = [a(:,1)./r a(:,2)./r a(:,3)./r];end
function r = kk_mx_vec_mult(a,b)
%%% a is nx9, b is nx3, output r is nx3
r = [ (a(:,1).*b(:,1)+ a(:,4).*b(:,2)+ (a(:,7).*b(:,3))) ...
      (a(:,2).*b(:,1)+ a(:,5).*b(:,2)+ (a(:,8).*b(:,3))) ...
      (a(:,3).*b(:,1)+ a(:,6).*b(:,2)+ (a(:,9).*b(:,3))) ];


































