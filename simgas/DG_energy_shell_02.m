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

%% calculate contravariant basis
m_u1 = kk_cross(Xp,n)./SSn(:,ones(1,3));
m_u2 = kk_cross(n,Xt)./SSn(:,ones(1,3));
% % %% calculate the contravariant metric tensor components g_uij
% % g_u11 = kk_dot(m_u1, m_u1);
% % g_u12 = kk_dot(m_u1, m_u2);
% % g_u22 = kk_dot(m_u2, m_u2);

%% test that g = dyadic(m_li, m_ui);  
m_l1 = Xt;
m_l2 = Xp;
m_l3 = n;

m_u3 = n;
% g_l11 = E;
% g_l22 = G;
% g_l12 = F;
% g_l13 = kk_dot(m_l1, m_l3);
% g_l23 = kk_dot(m_l2, m_l3);
% g_l33 = kk_dot(m_l3, m_l3);
% 
% g_u33 = kk_dot(m_u3, m_u3);
% g_u13 = kk_dot(m_u1, m_u3);
% g_u23 = kk_dot(m_u2, m_u3);


% m_u1_t = g_u11(:,ones(1,3)).*m_l1 + g_u12(:,ones(1,3)).*m_l2;   % compare to m_u1 --- passed
% m_u2_t = g_u12(:,ones(1,3)).*m_l1 + g_u22(:,ones(1,3)).*m_l2;   % compare to m_u1 --- passed

% % p = 100;
% % g_l = [[g_l11(p,:) g_l12(p,:) g_l13(p,:)];[g_l12(p,:) g_l22(p,:) g_l23(p,:)]; [g_l13(p,:) g_l23(p,:) g_l33(p,:)] ];
% % % disp(g_l);
% % g_u = [[g_u11(p,:) g_u12(p,:) g_u13(p,:)];[g_u12(p,:) g_u22(p,:) g_u23(p,:)]; [g_u13(p,:) g_u23(p,:) g_u33(p,:)] ];
% % % disp(g_u);
% % 
% % m_l1t = Xt(p,:);
% % m_l2t = Xp(p,:);
% % m_l3t = n(p,:);
% % 
% % m_u1t = Xt(p,:);
% % m_u2t = Xp(p,:);
% % m_u3t = n(p,:);
% % 
% % g_p = (m_l1t'* m_u1t) + (m_l2t'* m_u2t) + (m_l3t'* m_u3t);
% % % disp('gp');
% % % disp(g_p);

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

%% calculate shape energy
Eb = 1/2*sum((2.*H(:)-2.*C_o(:)).^2.*w.*SSn)./8/pi;   % reduced bending energy
% E = Eb;
if isempty(phys.SSn_o),     % i.e. in case we want to calculate the undeformed contributions
    E_stretch = 0;
    phys.Xt_o = Xt;
    phys.Xp_o = Xp;
    phys.E_o  = E;
    phys.G_o  = G;
    phys.F_o  = F;
    phys.L_o = L;
    phys.M_o = M;
    phys.N_o = N;
    phys.SS_o = SS;
    phys.SSn_o = SSn;
    phys.m_l1 = m_l1;
    phys.m_l2 = m_l2;
    phys.m_l3 = m_l3;
    phys.m_u1 = m_u1;
    phys.m_u2 = m_u2;
    phys.m_u3 = m_u3;
else                        % i.e. in case we really want the shear energy relative to the undeformed shape
    ai = (SSn./phys.SSn_o - 1);     % area_invariant
    %%% calculate stretch energy
    E_stretch = phys.k_stretch/2 * ...
        sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
    
    %%% Calculate  DGT = kron(defgi, undefgi);
    DGT = kk_kron(m_l1,phys.m_u1) + kk_kron(m_l2,phys.m_u2) + kk_kron(m_l3, phys.m_u3);
    FT  = kk_transpose(DGT);
    FTF = kk_mx_mult(FT,DGT);
    %%% calculate the material (Lagrangian) strain
    eps = 1/2*FTF;
    eps(:,1) = eps(:,1)-1/2;eps(:,5) = eps(:,5)-1/2;eps(:,9) = eps(:,9)-1/2;
    eps11 = eps(:,1);
    eps22 = eps(:,5);
    eps12 = eps(:,2);
    %%% put the tensors in their PAS
    b = -(eps11+eps22); c = eps11.*eps22 -eps12.*eps12; % a is always 1
    E1 = (-b + sqrt(b.*b - 4.* c))./2;E2 = (-b - sqrt(b.*b - 4.* c))./2;%disp([E1 E2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta  = (E1 + E2 - ai)./(1 + ai);%b1 = 0.7; b2 = 0.75;
    %Calculate shear energy
    E_shear = phys.k_shear *sum((beta + phys.b1 * ai .*beta + phys.b2 * beta.*beta).*phys.SSn_o.*w) ;
    %%disp(E_shear);
end

% % %%% calculate the shear energy based on triangles
% % X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
% % [E_shear] = shear_calc(X,phys.C_slave,phys.pre,phys.b1, phys.b2, phys.k_shear);
% % disp(E_shear);
E = full(Eb + E_shear + E_stretch);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the deformation gradient tensor DGT (assuming constant
%%%% thickness)
%%%% DGT ~=(g + x_3 * kappa) * (m_alpha (x) undeformed_m_alpha) * (undeformed_g - x_3 * undeformed_kappa)
%%%% DGT ~= shear + bending + thickness_change


%%% define the metric tensors directly from the coefficients of the first fundamental form
% % % % % % % zz = zeros(length(G),1);
% % % % % % % on = ones(length(G), 1);
% % % % % % % g_u = [phys.E_o phys.F_o zz phys.F_o phys.G_o zz zz zz on]; % for the undeformed configuration (passed through phys)
% % % % % % % g_d = [E F zz F G zz zz zz on]; % for the deformed configuration
% % % % % % % 
% % % % % % % 
% % % % % % % %%% define the metric tensors as the addition of three outer products
% % % % % % % g_u = kk_kron(phys.m_l1,phys.m_u1) + kk_kron(phys.m_l2,phys.m_u2) + kk_kron(phys.m_l3, phys.m_u3);
% % % % % % % g_d = kk_kron(m_l1,m_u1) + kk_kron(m_l2,m_u2) + kk_kron(m_l3, m_u3);
% % % % % % % 
% % % % % % % %%% 
% % % % % % % k_u = [(phys.L_o.*phys.G_o - phys.M_o.*phys.F_o) (phys.M_o.*phys.E_o -phys.L_o.*phys.F_o) zz (phys.M_o.*phys.G_o-phys.N_o.*phys.F_o) (phys.N_o.*phys.E_o-phys.M_o.*phys.F_o) zz zz zz on]./phys.SSn_o(:,ones(1,9));
% % % % % % % k_d = [(L.*G - M.*F) (M.*E -L.*F) zz (M.*G-N.*F) (N.*E-M.*F) zz zz zz on]./SSn(:,ones(1,9));
% % % % % % % % calculate the DGT terms; each term is an nx9 array
% % % % % % % DGT1 = [g_d + phys.D .*k_d];% the first term in the definition of the DGT 
% % % % % % % DGT2 = kk_kron(m_l1, m_u1) + kk_kron(m_l2, m_u2);
% % % % % % % DGT3 = [g_u - phys.D .*k_u];% approximation of the third term in the definition of the DGT 
% % % % % % % 
% % % % % % % DGT = kk_mx_mult(DGT2, DGT3);
% % % % % % % DGT = kk_mx_mult(DGT1, DGT);
% % % % % % % 
% % % % % % % detF = kk_mx_det(DGT);
% % % % % % % trFTF = kk_trace(kk_mx_mult(kk_transpose(DGT), DGT));
% % % % % % % 
% % % % % % % E = phys.lambda/2*log(detF.^2) + phys.miu/2.*(trFTF-3)-phys.miu/2.*log(detF);
% % % % % % % E = sum(E.*w.*SSn);


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
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      















