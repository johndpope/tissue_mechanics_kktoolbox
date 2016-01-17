function [E, A, V, v_red, AreaDiff, Eb, EADE, E_stretch, E_shear, T, H, phys] = ...
         DG_energy(X_o,phys)

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

cm_s  = [xclks(1) yclks(1) zclks(1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%% Calculate the geometrically important quantities
X  =[x(:) y(:) z(:)];
Xt = [xt(:) yt(:) zt(:)];
Xp = [xp(:) yp(:) zp(:)];
Xpp = [xpp(:) ypp(:) zpp(:)];
Xtp = [xtp(:) ytp(:) ztp(:)];
Xtt = [xtt(:) ytt(:) ztt(:)];

E = kk_dot(Xt,Xt);
F = kk_dot(Xt,Xp);
G = kk_dot(Xp,Xp);
SS = (kk_cross(Xt,Xp));SSn = sqrt(E.*G-F.*F);
n = SS./SSn(:,ones(1,3));
L = kk_dot(Xtt,n);
M = kk_dot(Xtp,n);
N = kk_dot(Xpp,n);

%% calculation of geometric properties for current configuration
V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
A = sum(w.*SSn);
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
K = (L.*N - M.*M)./(E.*G-F.*F);
h = 1./A.*sum(H(:).*w.*SSn);
T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
Eb = 1/2*sum((2.*H(:)-2.*C_o(:)).^2.*w.*SSn)./8/pi;   % reduced bending energy
Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;
v_red = V/Vo;

% E = Eb;
if isempty(phys.SSn_o),     % i.e. in case we want to calculate the undeformed contributions
    E_stretch = 0;
    phys.Xt_o = Xt;
    phys.Xp_o = Xp;
    phys.E_o  = E;
    phys.G_o  = G;
    phys.F_o  = F;
    phys.SS_o = SS;
    phys.SSn_o = SSn;
else                        % i.e. in case we really want the shear energy relative to the undeformed shape
    %%% area invairiant using differential area element
    ai = (SSn./phys.SSn_o - 1); % area invariant
    %%% calculate stretch energy
    E_stretch = phys.k_stretch/2 * ...
        sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
end

%%% calculate the shear energy based on triangles
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
[E_shear] = shear_calc(X,phys.C_slave,phys.pre,phys.b1, phys.b2, phys.k_shear);

E = full(Eb + E_shear + E_stretch);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















