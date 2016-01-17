function [E, A, V, v_red, AreaDiff, Eb, EADE, E_stretch, E_shear, T, H, phys] = ...
    DG_energy_shell(X_o,phys)
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
%cm_m   = phys.cm;

%k_stretch= phys.k_stretch;

% %%%% obtain quantities for undeformed configuration
% E_u = phys.E_u;
% F_u = phys.F_u;
% G_u = phys.G_u;
% L_u = phys.L_u;
% M_u = phys.M_u;
% N_u = phys.N_u;
% SSn_u = phys.SSn_u;

X_o = recover_parms(X_o,phys);
X_o = over_Basis(X_o);
[xclks yclks zclks] = get_xyz_clks(X_o);
%cm_s  = [xclks(1) yclks(1) zclks(1)];
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
m_l3 = n;

bb = kk_dot(m_l1, kk_cross(m_l2,m_l3));bb = bb(:);  % formulation as in AMoS p.658
bb = kk_norm(kk_cross(m_l1,m_l2));                  % bb is essentially SSn
m_u1 = kk_cross(Xp,n)./bb(:,ones(1,3));
m_u2 = kk_cross(n,Xt)./bb(:,ones(1,3));
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
g_u33 = kk_dot(m_u3, m_u3);
g_u13 = kk_dot(m_u1, m_u3);
g_u23 = kk_dot(m_u2, m_u3);
%% calculate the coefficients of the second fundamental form
L = kk_dot(Xtt,n);
M = kk_dot(Xtp,n);
N = kk_dot(Xpp,n);
%% calculate  geometric properties: Area, Volume, reduced volume, Curvatures and self-check total Gaussian curvature
V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
A = sum(w.*SSn);
Vo = 4/3*pi*(A/4/pi)^(3/2);
v_red = V/Vo;
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
K = (L.*N - M.*M)./(E.*G-F.*F);
h = 1./A.*sum(H(:).*w.*SSn);
T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
%% calculate bending energy
Eb = 1/2*sum((2.*H(:)-2.*C_o(:)).^2.*w.*SSn)./8/pi;   % reduced bending energy
%% calculate in-plane shear strain as in AoPaS page 214 ---- needs work
omega = kk_dot(tt, tp);
b5 = omega;

%% Works!! Calculate midplane lagrange strain tensor as in http://www.scribd.com/doc/27473029/5/Strain-in-Orthogonal-Curvilinear-Coordinates
%%%%     originally from AMoS page 667 but didn't work: still figuring out why this is correct (normalization skipped in AMoS?)
%%% using loops!!
lambda = zeros(length(Xt),2);
gamma11 = 1/2*(g_l11-phys.g_l11)./sqrt(phys.g_l11)./sqrt(phys.g_l11);
gamma12 = 1/2*(g_l12-phys.g_l12)./sqrt(phys.g_l11)./sqrt(phys.g_l22);
gamma21 = 1/2*(g_l12-phys.g_l12)./sqrt(phys.g_l11)./sqrt(phys.g_l22);
gamma22 = 1/2*(g_l22-phys.g_l22)./sqrt(phys.g_l22)./sqrt(phys.g_l22);
gamma13 = 1/2*(g_l13-phys.g_l13)./sqrt(phys.g_l11)./sqrt(phys.g_l33);
gamma23 = 1/2*(g_l23-phys.g_l23)./sqrt(phys.g_l22)./sqrt(phys.g_l33);
gamma31 = 1/2*(g_l13-phys.g_l13)./sqrt(phys.g_l33)./sqrt(phys.g_l11);
gamma32 = 1/2*(g_l23-phys.g_l23)./sqrt(phys.g_l33)./sqrt(phys.g_l22);
gamma33 = 1/2*(g_l33-phys.g_l33)./sqrt(phys.g_l33)./sqrt(phys.g_l33);
eps_vec = [gamma11 gamma12 gamma13 gamma21 gamma22 gamma23 gamma31 gamma32 gamma33];    % construct the strain tensor
for ix = 1:length(Xt)
    eps = reshape(eps_vec(ix,:), 3,3);
    FTF = 2*eps+eye(3,3);             % Right Cauchy-Green tensor
    eFTF= eig(FTF);                   % calculate eigenvalues of right Cauchy-Green tensor (square root of which provides the principal stretches AMoS p.29)
    lambda(ix,:) = [sqrt(eFTF(1)) sqrt(eFTF(2))];       % store the principal stretches
end
b4 = (lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2); % shear invariant
ai = lambda(:,1).*lambda(:,2) - 1;     % area dilation invariant (as function of principal stretches
E_stretch4 = phys.k_stretch/2 * ...
    sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% area invariant DG
ai = (SSn./phys.SSn_o - 1);     % area_invariant
E_stretch1 = phys.k_stretch/2 * ...
    sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area

%% % calculate the shear and stretch energy based on triangles
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
[E_shear4 E_stretch3 E1 E2 b_tri] = shear_calc(X,phys.C_slave,phys.pre,phys.a3, phys.a4,phys.b1, phys.b2, phys.k_shear);
%% % calculate Strain tensor components from first term of F in AMoS: p 665 
%%%  using loops
lambda = zeros(length(Xt),2);
beta = zeros(length(Xt),1);
E1v = zeros(length(Xt),1);
E2v = zeros(length(Xt),1);
eigFTF = zeros(length(Xt), 3);
DGT_vec = zeros(length(Xt), 9);
for ix = 1:length(Xt)
    DGT = kk_kron(m_l1(ix,:), phys.m_u1(ix,:)) +...
          kk_kron(m_l2(ix,:), phys.m_u2(ix,:)) +...
          kk_kron(m_l3(ix,:), phys.m_u3(ix,:));% +...
%           kron(m_l1(ix,:), phys.m_u2(ix,:)) +...
%           kron(m_l2(ix,:), phys.m_u1(ix,:)) + ...
%           kron(m_l1(ix,:), phys.m_u3(ix,:)) + ...
%           kron(m_l3(ix,:), phys.m_u1(ix,:)) + ...
%           kron(m_l2(ix,:), phys.m_u3(ix,:)) + ...
%           kron(m_l3(ix,:), phys.m_u2(ix,:));
    DGT_vec(ix,:) = DGT(:)';
    DGT = reshape(DGT,3,3);%     DGT = DGT(1:2,1:2);
    FTF = DGT'*DGT;     % right Cauchy-Green tensor
    eFTF= eig(FTF);
    lambda(ix,:) = [sqrt(eFTF(1)) sqrt(eFTF(2))];
    eigFTF(ix,:) = eFTF(:)';
    eps = 1/2*(FTF-eye(size(FTF)));
    Eig_vec = eig(eps, 'nobalance');
    E1 = Eig_vec(1);
    E2 = Eig_vec(2);
    beta(ix) = (E1 + E2 - ai(ix))./(1 + ai(ix));    % shear component
    E1v(ix) = E1;
    E2v(ix) = E2;
end
b1 = beta;
b2 = (lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2); % expression for beta as function of principal stretches

% stretch 2
ai = lambda(:,1).*lambda(:,2) - 1;     % area_invariant (as function of principal stretches
E_stretch2 = phys.k_stretch/2 * ...
    sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
%% % calculate strain tensor vectorized version
    %%% Calculate  DGT = kron(defgi, undefgi);
    DGT = kk_kron(m_l1,phys.m_u1) + kk_kron(m_l2,phys.m_u2) + kk_kron(m_l3, phys.m_u3);
    FT  = kk_transpose(DGT);
    C = kk_mx_mult(FT,DGT);
    FTF = C;        % right Cauchy-Green deformation tensor
    
    %%% calculate the material (Lagrangian) strain
    eps = 1/2*FTF;
    eps(:,1) = eps(:,1)-1/2;eps(:,5) = eps(:,5)-1/2;eps(:,9) = eps(:,9)-1/2;
    eps11 = eps(:,1);
    eps22 = eps(:,5);
    eps12 = eps(:,2);
%% % put the tensors in their PAS
b = -(eps11+eps22); c = eps11.*eps22 -eps12.*eps12; % a is always 1
E1 = (-b + sqrt(b.*b - 4.* c))./2;E2 = (-b - sqrt(b.*b - 4.* c))./2;%disp([E1 E2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ai = (SSn./phys.SSn_o - 1);     % area_invariant
b3  = (E1 + E2 - ai)./(1 + ai);%b1 = 0.7; b2 = 0.75;
%%
%%%Calculate shear energy
ai = (SSn./phys.SSn_o - 1);     % area_invariant
beta = b1;E_shear1 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
beta = b2;E_shear2 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
beta = b3;E_shear3 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
beta = b4;E_shear5 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
beta = b5;E_shear6 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;

%% reporting
disp('---------Stretch energy------------');
disp(['Based on original area invariant DG: ' num2str(E_stretch1)]);
disp(['AMoS p665 first term of F (using loops) using principal stretches: ' num2str(E_stretch2)]);
disp(['Based on triangles: ' num2str(E_stretch3)]);
disp(['Based on midplane lagrange strain tensor (normalized expression): ' num2str(E_stretch4)]);

disp('-------- Shear energy --------------');
disp(['AMoS p665 first term of F (using loops) -- beta calculated from strain tensor with ai from DG:  ' num2str(E_shear1)]);
disp(['AMoS p665 first term of F (using loops) -- beta calculated from principal stretches directly:  ' num2str(E_shear2)]);
disp(['AMoS p665 first term of F (vectorized) -- beta calculated from strain tensor with ai from DG:  ' num2str(E_shear3)]);
disp(['Based on triangles:  ' num2str(E_shear4)]);
disp(['Based on midplane lagrange strain tensor (normalized expression): ' num2str(E_shear5)]);
disp(['in-plane shear strain as in AoPaS page 214: ' num2str(E_shear6)]);
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
a1 = a(:,1);
b1 = a(:,2);
c1 = a(:,3);
d1 = a(:,4);
e1 = a(:,5);
f1 = a(:,6);
g1 = a(:,7);
h1 = a(:,8);
i1 = a(:,9);

a2 = b(:,1);
b2 = b(:,2);
c2 = b(:,3);
d2 = b(:,4);
e2 = b(:,5);
f2 = b(:,6);
g2 = b(:,7);
h2 = b(:,8);
i2 = b(:,9);


r = [ a1.*a2 + b1.*d2 + c1.*g2 ...
    a1.*b2 + b1.*e2 + c1.*h2 ...
    a1.*c2 + b1.*f2 + c1.*i2 ...
    a2.*d1 + d2.*e1 + f1.*g2 ...
    b2.*d1 + e1.*e2 + f1.*h2 ...
    c2.*d1 + e1.*f2 + f1.*i2 ...
    a2.*g1 + d2.*h1 + g2.*i1 ...
    b2.*g1 + e2.*h1 + h2.*i1 ...
    c2.*g1 + f2.*h1 + i1.*i2];
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
function r = kk_norm(a)
%%%% a is a nx3 vector
%%%% r is a nx1 vector of the norm of every row in a
r = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);




































