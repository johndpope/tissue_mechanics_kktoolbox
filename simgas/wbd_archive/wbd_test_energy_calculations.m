clc;
% load X_o_dros_embryo_04.mat;
s = shp_surface;X_o = s.X_o;
L_max = 5;
gdim = 60;
model_fun = @DG_energy_shell_09;
constraints_fun = @DG_constraints;
%%%% prepare slave surface
s_slave = shp_surface(L_max,gdim);
s_slave.X_o = X_o;
s_slave = update(s_slave);
As = s_slave.A;Vs = s_slave.V;
Xs = X_o;           % slave surface
X_o = shp_surface.tr(X_o, L_max);
X_o = times_basis(X_o);
%% prepare physical parameters
phys.D  = 15 * 1e-6;  % units micron
phys.A_o = As;
phys.V_o = Vs;
phys.x3 = 1;
% phys.Young = 1000;
% phys.Poiss = 0.3;

phys.lambda = 577;   % N m^-2
phys.miu = 384;      % N m^-2


phys.kb = 1;%2e-19; %units in joule
phys.k_stretch = 1;%1e-2;%1e-20;%   set to artificially weak membrane to keep the triangles from shearing. was 5e-18;
phys.k_shear = 1;%phys.k_stretch/2;
phys.a3 = -2; 
phys.a4 = 8;
phys.b1 = 0.7; 
phys.b2 = 0.75;
%% prepare triangle-based energy calculation
nicos = 4;      % number of icosahedron subdivisions used for self-intersection test
[X, phys.C_slave]=surface_mesh.sphere_mesh_gen(nicos);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
%%% generate the self-intersection test basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
phys.Y_LK_self  = zeros(N, M, 'single');
for S = 1:length(L),
    phys.Y_LK_self(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
[xclks yclks zclks] = get_xyz_clks(over_Basis(X_o));
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
%%% precalculation for using shear with triangles
phys.pre = shear_stretch_precalc(X,phys.C_slave);    % precalcualtion for fast vectorized calculation of shear using triangles
%% prepare the DG energy calculation
[X_temp B Y_ge] = DG_precalc(phys.A_o, L_max,L_max, gdim);phys.B = B;

phys.C_o = 0;
phys.SSn_o = [];
phys.uX_o = X_o;    % undeformed configuration -- to get the correct one use overbasis
[fit_parms phys]= parmvec_gen(X_o,phys);
[phys] = DG_energy_shell_prep(fit_parms,phys);
%[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H, phys] = model_fun(fit_parms,phys);   % test the energy calculation
%disp(E);
%%
% % %% % change the shape
% s = shp_surface('bowling_pin');
%  X_o = tr(s.X_o, L_max);
X_o = X_o*10;
X_o(7) = 5;

[fit_parms phys]= parmvec_gen(X_o,phys);
[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);   % test the energy calculation

% disp(E);




















