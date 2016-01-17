function [area_percent] = check_MP_area(X_o,phys)
%% process input
Y_P = phys.B.Y_P;
Y_T = phys.B.Y_T;
p    = phys.B.p;
t    = phys.B.t;
w    = phys.B.w;
X_o = recover_parms(X_o,phys);
X_o = over_Basis(X_o);
[xclks yclks zclks] = get_xyz_clks(X_o);
%% generate x y z, and their individual derivatives with respect to theta and phi
gdimp = length(p);gdimt = length(t);
c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
xp =(sum(c.*Y_P,3));
xt = (sum(c.*Y_T,3));

c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
yp =(sum(c.*Y_P,3));
yt = (sum(c.*Y_T,3));

c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
zp =(sum(c.*Y_P,3));
zt = (sum(c.*Y_T,3));
%% generate X surface vector and derivatives of X
Xt = [xt(:) yt(:) zt(:)];       % m_l1
Xp = [xp(:) yp(:) zp(:)];       % m_l2
%% calculate coefficients of the 1st fundamental form (components of g_lij)
E = kk_dot(Xt,Xt);          % g_l11
F = kk_dot(Xt,Xp);          % g_l12 and g_l21
G = kk_dot(Xp,Xp);          % g_l22
%% calculate differential area element and surface normal

SSn = sqrt(E.*G-F.*F);
%% calculate  geometric properties: Area, Volume, reduced volume, Curvatures and self-check total Gaussian curvature
A = sum(w.*SSn);

A_MP = sum(w.*SSn.*phys.MP);          % area covered by mesoderm primordium
area_percent = A_MP/A * 100;
disp(['MP area percent: ' num2str(area_percent)]);












