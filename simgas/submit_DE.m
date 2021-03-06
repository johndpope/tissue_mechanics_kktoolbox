function [X_o_res_vec, E_vec, exitflag_vec] = submit_DE(in)
%%%% gets called from a wbd_meta_simgas_XXX script or from a GUI
%%% [1] pre-process the input required for gastrulation simulations
%%% [2] submit the calculation to the optimizer for every calculation
%%%     sequentially
%%% [3] Save the results
%%% [4] plot the results and save high res movies/images if required
%% interpret input struct "in"
L_max = in.L_max;           % maximum order of the spherical harmonics parameterization
L_max_gef = in.L_max_gef;
gdim = in.gdim;
optim_method = in.optim_method;
nicos = in.nicos;
nicovm = in.nicovm;
pos = in.pos;
model_fun = in.model_fun;
constraints_fun = in.constraints_fun;
phys = in.phys;

%%% calculate  material properties (given Young's modulus and Poisson ratio)
phys.lambda = phys.Poiss*phys.Young/(1-2*phys.Poiss)/(1+phys.Poiss); % Lame's first parameter% phys.lambda = 577;   % N m^-2
phys.miu = phys.Young/2/(1+phys.Poiss);  % Shear modulus, Lame's second parameter (also G sometimes)% phys.miu = 385;      % N m^-2
phys.bulk = phys.Young/3/(1-2*phys.Poiss);
phys.kb = phys.Young * phys.D^3/12/(1-phys.Poiss^2); %
% phys.k_shear = phys.miu;
% phys.k_stretch = phys.k_shear * 2;

%%% Material parameters for mesoderm primordium MP
phys.YoungMP    = phys.Young*phys.MPfac;
phys.PoissMP    = phys.Poiss;
phys.lambdaMP   = phys.PoissMP*phys.YoungMP/(1-2*phys.PoissMP)/(1+phys.PoissMP);
phys.miuMP      = phys.YoungMP/2/(1+phys.PoissMP);
phys.bulkMP     = phys.YoungMP/3/(1-2*phys.PoissMP);
phys.kbMP       = phys.MPfac*phys.kb;

% % phys.del_indx = [];
if isempty(in.fn_del_ind), phys.del_indx = [];
else load(in.fn_del_ind, 'del_indx');phys.del_indx = del_indx;
end
%% initialize morphologies (master, slave, starting shape)
load(in.fn_master,'X_o');X_o = shp_surface.tr(X_o, L_max);
phys.X_o_vm = tr(shp_scale(X_o, in.VMscale), L_max);  % the geometry of the vitteline membrane
s_master = shp_surface(L_max,gdim);
s_master.X_o = phys.X_o_vm;
s_master = update(s_master);
s_master.use_camorbit = 0;
%%%% prepare slave surface --- i.e. the embryo tissue
load(in.fn_slave,'X_o');X_o = shp_surface.tr(X_o, L_max);
s_slave = shp_surface(L_max,gdim);
s_slave.X_o = X_o;
s_slave = update(s_slave);
Xs = X_o;           % slave surface
As = s_slave.A;
Vs = s_slave.V;
phys.A_o = As;
phys.V_o = Vs;
if isfield(phys,'X_o_reference'),
    disp('Found (and using) input reference configuration!');
    X_o_reference = times_basis(phys.X_o_reference);
else
    disp('using slave configuration as reference shape');
    X_o_reference = times_basis(X_o);
end
%%%% prepare starting surface 
load(in.fn_start,'X_o');X_o = shp_surface.tr(X_o, L_max);
X_o = times_basis(X_o); % starting shape
%% initialize basis functions
[X_temp B Y_ge] = DG_precalc(phys.A_o, L_max,L_max_gef, gdim);phys.B = B; 
gdimp = length(B.p);gdimt = length(B.t);
%% initialize gene expression fields
load(in.fn_ge,'gtw', 'gsnl', 'ghkb');
gtw = sh_surface.tr_xc(gtw,L_max_gef);
gsnl = sh_surface.tr_xc(gsnl,L_max_gef);
ghkb = sh_surface.tr_xc(ghkb,L_max_gef);
%%% in case of ge.mat
% s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s.xc = gtw;s = sh_rot(s,0, 0, pi);gtw = s.xc;
% s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s.xc = gsnl;s = sh_rot(s,0, 0, -pi/12);gsnl = s.xc;

if(in.rotation_correction),
    %%% in case of ge_latest.mat we need to perform some rotation corrections
    s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s.xc = gtw;s = sh_rot(s,0, 0, pi/2);gtw = s.xc;
    s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s.xc = gsnl;s = sh_rot(s,0, 0, pi/2);gsnl = s.xc;
    s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s.xc = ghkb;s = sh_rot(s,0, 0, pi/2);ghkb = s.xc;   
end
c = gtw(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);gtwf = sum(c.*Y_ge,3);
c = gsnl(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);gsnlf = sum(c.*Y_ge,3);
c = ghkb(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);ghkbf = sum(c.*Y_ge,3);
%% initialize self-intersection test
[X, phys.C_slave]=surface_mesh.sphere_mesh_gen(nicos);
load(in.fn_data_self,'TP');
phys.TP_self = TP;
[res_self] = tri_tri_self_intersect(X,phys.C_slave, uint16(phys.TP_self));

[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
%%% generate the self-intersection test basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
phys.Y_LK_self  = zeros(N, M, 'single');
for S = 1:length(L),
    phys.Y_LK_self(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
%% initialize master-intersection test
[X2, phys.C_vm]=surface_mesh.sphere_mesh_gen(nicovm);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
%%% generate the vm basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
Y_LK_vm  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK_vm(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
[xc yc zc] = get_xyz_clks(phys.X_o_vm);
phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for vitteline membrane
phys.cm = sum(phys.X_vm,1)./length(phys.X_vm);
load(in.fn_tdb_data,'TP');phys.TP_vm = TP;
[res] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
%% calculate initial energy to define reference shape properties
% [xclks yclks zclks] = get_xyz_clks(over_Basis(X_o));
% X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
phys.MP = 0;
phys.miu_eff = 0;
phys.bulk_eff = 0;
phys.kb_eff = 0;
phys.k_shear_eff = 0;
phys.k_stretch_eff = 0;
phys.C_o = 0;
phys.uX_o = X_o;    % undeformed configuration -- to get the correct one use overbasis
[fit_parms phys]= parmvec_gen(X_o_reference,phys);
[phys] = DG_energy_shell_prep(fit_parms,phys);  % for the calculation of the undeformed shape geometrical parameters
[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);   % test the energy calculation
phys.C_o = H;       % to make sure that the shape has zero curvature energy if it were the reference shape
%% prepare the effect of the gene expression field
gsnlf = gsnlf/norm(gsnlf(:));
gtwf = gtwf/norm(gtwf(:));
ghkbf = ghkbf/norm(ghkbf(:));
facf = gsnlf(:).^(phys.psnl).*gtwf(:).^(phys.ptw) .*ghkbf(:).^(phys.phkb);     % define how individual genes affect the curvature
facf = facf./norm(facf(:));
phys.C_o = H + phys.lambda_ge.*facf;      % modulate the preferred curvature
% phys.C_o = H + phys.lambda_ge;
disp(['Using curvature ' num2str(phys.lambda_ge) ' for MP region']);
%%%% define MP region
phys.MP = facf;
phys.MP(phys.MP<phys.MP_cutoff*max(phys.MP)) = 0;
phys.MP(phys.MP>0) = 1; % defines MP region
phys.MP = logical(phys.MP);
% phys.C_o = H + phys.lambda_ge*phys.MP;
%%%% generate the material properties that will be used for each region
miu = phys.miu.*ones(size(facf));miu(phys.MP) = phys.miuMP;
bulk = phys.bulk.*ones(size(facf));bulk(phys.MP) = phys.bulkMP;
kb = phys.kb.*ones(size(facf));kb(phys.MP) = phys.kbMP;
% k_shear = phys.k_shear.*ones(size(facf));
% k_shear(phys.MP) = phys.k_shearMP;
% k_stretch = phys.k_stretch.*ones(size(facf));
% k_stretch(phys.MP) = phys.k_stretchMP;

phys.miu_eff = miu;
phys.bulk_eff = bulk;
phys.kb_eff = kb;
% phys.k_shear_eff = k_shear;
% phys.k_stretch_eff = k_stretch;
%% plot MP or facf
%hist([gsnlf(:)/norm(gsnlf(:)) gtwf(:)/norm(gtwf(:)) ghkbf(:)/norm(ghkbf(:))])
% plot_MP(fit_parms, double(phys.MP),phys);
subplot(5,1,1);plot_MP_2d(fit_parms, double(gsnlf./max(gsnlf(:))),phys, 'snail');
subplot(5,1,2);plot_MP_2d(fit_parms, double(gtwf./max(gtwf(:))),phys, 'twist');
subplot(5,1,3);plot_MP_2d(fit_parms, double(ghkbf./max(ghkbf(:))),phys, 'hkb');
subplot(5,1,4);plot_MP_2d(fit_parms, double(facf./max(facf(:))),phys, 'effective gene activity');
subplot(5,1,5);plot_MP_2d(fit_parms, double(phys.MP),phys, 'MP region');

%% plot perspective views
plot_MP(fit_parms, double(gsnlf./max(gsnlf(:))),phys);axis off;axis equal
view(-90,0);lighting flat;camlight;
view(-270,-90);camlight;axis off
% plot_MP(fit_parms, double(gtwf./max(gtwf(:))),phys, 'twist');
% plot_MP(fit_parms, double(ghkbf./max(ghkbf(:))),phys, 'hkb');
% plot_MP(fit_parms, double(facf./max(facf(:))),phys, 'effective gene activity');
% plot_MP(fit_parms, double(phys.MP),phys, 'MP region');

%% plot gene expression coefficients histogram
figure;
subplot(3,1,1);hist(gsnl,100);title('snail');xlim([-60 60]);ylim([0 180]);
subplot(3,1,2);hist(gtw,100);title('twist');xlim([-60 60]);ylim([0 180]);
subplot(3,1,3);hist(ghkb,100);title('huckebein');xlim([-60 60]);ylim([0 180]);
%%
phys.V_o = V * 1.00;       % define the volume constraint
phys.A_o = A;
Vo = 4/3*pi*(A/4/pi)^(3/2);
phys.v_red_o = V/Vo;
fit_parms = parmvec_gen(X_o,phys);
[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, H] = model_fun(fit_parms,phys);   % test the energy calculation
[c ceq] = constraints_fun(fit_parms,phys);   % test the energy calculation
str = sprintf('Reduced volume: \t%.4f\n Area: \t\t\t%.4f\n Volume: \t%.4f\nTotal E.:\t%.4f\nBending E.: \t%.4f\nADE E.: \t%.4f\n Stretch E.: %.4f\n Shear E.: %.4f',...
    V./(4/3*pi*(A/4/pi).^(3/2)), A, V, E, Eb, EADE, E_stretch, E_shear);disp(str);
[mp_area] = check_MP_area(fit_parms,phys);
%% optimize
[X_o_res_vec, E_vec, exitflag_vec] = combi_optim_02(X_o,phys,s_master,pos, Y_LK_vm, in.mcmc_options,in.newton_options, model_fun, constraints_fun);
%% save and plot output
str = sprintf('save %s.mat in L_max gdim X_o_res_vec X_o_reference s_master phys pos Y_LK_vm exitflag_vec;',in.runName);eval(str);
% plot_time_sequence_pretty_03(X_o_res_vec(1:end,:),phys,[]);
% plot_time_sequence_pretty(X_o_res_vec(1:end,:),s_master,phys,pos, Y_LK_vm);