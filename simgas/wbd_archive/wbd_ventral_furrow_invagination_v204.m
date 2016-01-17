%%%% modeling ventral furrow invagination

clc;close all;global counter
counter = 0;
runName = 'DG_dros_embryo';
%% Step 0: configure numerical method for shape and gef representation
L_max = 7;           % maximum order of the spherical harmonics parameterization
L_max_gef = 5;
gdim = 30;
optim_method = 2;
pos = [ 0 0 0];
%% step 1: load embryo and configure master surface
%%%% prepare master surface
load X_o_dros_embryo_04.mat;X_o = shp_surface.tr(X_o, L_max);
phys.X_o_vm = tr(shp_scale(X_o, 1.05), L_max);  % the geometry of the vitteline membrane
s_master = shp_surface(L_max,gdim);
s_master.X_o = X_o;
s_master = update(s_master);
s_master.use_camorbit = 0;
Xm = s_master.X_o;    %
Am = s_master.A;
Vm = s_master.V;

%%%% prepare slave surface
s_slave = shp_surface(L_max,gdim);
s_slave.X_o = X_o;
s_slave = update(s_slave);
As = s_slave.A;Vs = s_slave.V;
Xs = X_o;           % slave surface


X_o = times_basis(X_o);

%% % L = 32;
% % % nico = 6;
% % % 
% % % s = sh_surface(L,sh_basis(L,120));
% % % s.xc = sh_surface.tr_xc(gtw,L);s = sh_rot(s,0, 0, pi);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'twist_globe.fig');myaa;print -dtiff -r600 twist_globe.tif;
% % % d = add_sf(d,'twist',s);
% % % 
% % % s.xc = sh_surface.tr_xc(gsnl,L);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'snail_globe.fig');myaa;print -dtiff -r600 snail_globe.tif;
% % % d = add_sf(d,'snail',s);
% % % 
% % % s.xc = sh_surface.tr_xc(ghkb,L);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'hkb_globe.fig');myaa;print -dtiff -r600 hkb_globe.tif;
% % % d = add_sf(d,'hkb',s);
% % % 
% % % plot_field(d,nico,1:3);view(24,-60);saveas(gcf,'overlayed.fig');myaa;print -dtiff -r600 overlayed.tif;
% % % plot_field(d,nico,1,1);view(24,-60);saveas(gcf,'twist.fig');myaa;print -dtiff -r600 twist.tif;
% % % plot_field(d,nico,2,2);view(24,-60);saveas(gcf,'snail.fig');myaa;print -dtiff -r600 snail.tif;
% % % plot_field(d,nico,3,3);view(24,-60);saveas(gcf,'hkb.fig');myaa;print -dtiff -r600 hkb.tif;
% % % save d_o d

%% Step 1: configure geometrical and physical quantities %%%%%
phys.D  = 15;  % units micron
phys.A_o = As;
phys.V_o = Vs;
phys.Young = 1000;
phys.Poiss = 0.3;


phys.kb = 2e-19; %units in joule
phys.k_stretch = 1e2;%1e-20;%   set to artificially weak membrane to keep the triangles from shearing. was 5e-18;
phys.k_shear = phys.k_stretch/2;
phys.a3 = -2; 
phys.a4 = 8;
phys.b1 = 0.7; 
phys.b2 = 0.75;



model_fun = @DG_energy;
constraints_fun = @DG_constraints;
%% Step 3: Precalculate quantities and test the energy calculation
[X_temp B Y_ge] = DG_precalc(phys.A_o, L_max,L_max_gef, gdim);phys.B = B; 
gdimp = length(B.p);gdimt = length(B.t);
%%% gene expression fields (snl, tw, hkb)
load ge;
gtw = sh_surface.tr_xc(gtw,L_max_gef);
gsnl = sh_surface.tr_xc(gsnl,L_max_gef);
ghkb = sh_surface.tr_xc(ghkb,L_max_gef);
% we need to rotate twist
s = sh_surface(L_max_gef,sh_basis(L_max_gef,gdim));s = sh_rot(s,0, 0, pi);gtw = s.xc;
c = gtw(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);gtw = sum(c.*Y_ge,3);
c = gsnl(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);gsnl = sum(c.*Y_ge,3);
c = ghkb(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);ghkb = sum(c.*Y_ge,3);

%% prepare quantities for self-intersection test
nicos = 2;      % number of icosahedron subdivisions used for self-intersection test
[X, phys.C_slave]=surface_mesh.sphere_mesh_gen(nicos);
%[res, phys.TP] = tri_tri_self_intersect(X,phys.C_self);  % generate the TP matrix for the fast self-intersection test
load DATA_self_nico_2.mat;  % loads precalculated TP matrix (located in the shp_surface class flolder)
phys.TP_self = TP;
[res_self] = tri_tri_self_intersect(X,phys.C_slave, uint16(phys.TP_self));

[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
%%% generate the self-intersection test basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
phys.Y_LK_self  = zeros(N, M, 'single');
for S = 1:length(L),
    phys.Y_LK_self(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
%% prepare quantities for vitteline membrane/embryo intersection test
nicovm = 2;      % number of icosahedron subdivisions used for vm for vm/embryo intersection test
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
load TDB_nico2.mat;phys.TP_vm = TP;
[res] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
%% % plot both vitteline membrane and embryo together
%[X1, C1, X2, C2] = shp_slice_x(over_Basis(X_o), 100);

[xclks yclks zclks] = get_xyz_clks(over_Basis(X_o));
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%

%%% precalculation for using shear with triangles
phys.pre = shear_stretch_precalc(X,phys.C_slave);    % precalcualtion for fast vectorized calculation of shear using triangles


%plot_shell(X, phys.C_slave,phys.D);
% patch('vertices',X,'faces',phys.C_slave,'FaceColor','red', 'FaceAlpha', 0.8);
% hold on;
% patch('vertices',phys.X_vm,'faces',phys.C_vm,'FaceColor','green', 'FaceAlpha', 0.6);
% view(3);axis equal;axis on;lighting phong;camlight;xlim([-300 300]);
%% test intersection calculations
[res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);% self-intersection test
[res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
%% test energy calculation
phys.C_o = 0;
phys.uX_o = X_o;    % undeformed configuration -- to get the correct one use overbasis
[fit_parms phys]= parmvec_gen(X_o,phys);
[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);   % test the energy calculation
phys.C_o = H;       % to make sure that the shape has zero energy as it is now
%% prepare gene expression distribution balance and test constraints calculation
% phys.C_o = H + gsnl(:)/max(gsnl(:));
lambda = 6e1;
phys.C_o = H+lambda*gsnl(:)./max(gsnl(:));
phys.V_o = V;       % fix the volume constraint to what it is at the beginning
phys.A_o = A;
Vo = 4/3*pi*(A/4/pi)^(3/2);
phys.v_red_o = V/Vo;
fit_parms = parmvec_gen(X_o,phys);
[E, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, H] = model_fun(fit_parms,phys);   % test the energy calculation
[c ceq] = constraints_fun(fit_parms,phys);   % test the energy calculation
str = sprintf('Reduced volume: \t%.4f\n Area: \t\t\t%.4f\n Volume: \t%.4f\nTotal E.:\t%.4f\nBending E.: \t%.4f\nADE E.: \t%.4f\n Stretch E.: %.4f\n Shear E.: %.4f',...
    V./(4/3*pi*(A/4/pi).^(3/2)), A, V, E, Eb, EADE, E_stretch, E_shear);disp(str);
%% Step 4: optimization
if optim_method == 0    % use *** Metropolis-Hastings *** naiive version
    %%%% configure MCMC method
    mcmc_options.TolCon          = 1e-3;
    mcmc_options.maxiter         = 1000;
    mcmc_options.sigQcon         = 1e-1;       % when constraints are not satisfied:standard deviation of Gaussian representing probability density for drawing samples
    mcmc_options.sigQsat         = 1e-3;       % case when constraints are satisfied
    mcmc_options.sigE            = 1.0e-3;       % standard deviation of the energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of MCMC configuration
    [X_o_res_vec, E_vec] = mcmc_optim(X_o,phys,s_master,pos, Y_LK_vm, mcmc_options, model_fun, constraints_fun);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif optim_method == 1
    options =   optimset(...
        'plotfcn',@DG_opt_plot,...
        'UseParallel','never',...
        'MaxFunEvals', 1e6,...
        'MaxIter', 1e4,...
        'DiffMaxChange', 1e1,...
        'DiffMinChange', 1e-2,...
        'DerivativeCheck','off',...
        'Algorithm','interior-point',...
        'AlwaysHonorConstraints','bounds',...
        'Hessian','lbfgs',...
        'InitTrustRegionRadius', 1,...
        'FinDiffType', 'central',...
        'ScaleProblem', 'none',...
        'GradObj','off',...
        'GradConstr','off',...
        'Display', 'iter',...
        'Diagnostics', 'off',...
        'FunValCheck', 'off',...
        'TolCon', 1e-3,...
        'TolFun', 1e-5,...
        'TolX', 1e-10,...
        'MaxSQPIter',1000);
    [X_o_res_vec, E_vec] = newton_optim(X_o,phys,s_master,pos, Y_LK_vm, options, model_fun, constraints_fun);
%%
elseif optim_method == 2
     %%% use a combination of MCMC and Newton methods
         %%%% configure MCMC method
    mcmc_options.TolCon          = 1e-3;
    mcmc_options.maxiter         = 5e4;
    mcmc_options.sigQcon         = 1e-2;       % when constraints are not satisfied:standard deviation of Gaussian representing probability density for drawing samples
    mcmc_options.sigQsat         = 1e-4;       % case when constraints are satisfied
    mcmc_options.sigE            = 1.0e-3;       % standard deviation of the energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of MCMC configuration
        options =   optimset('plotfcn',@DG_opt_plot,...
        'UseParallel','never',...
        'MaxFunEvals', 1e6,...
        'MaxIter', 1e4,...
        'DiffMaxChange', 1e1,...
        'DiffMinChange', 1e-2,...
        'DerivativeCheck','off',...
        'Algorithm','interior-point',...
        'AlwaysHonorConstraints','bounds',...
        'Hessian','lbfgs',...
        'InitTrustRegionRadius', 1,...
        'FinDiffType', 'central',...
        'ScaleProblem', 'none',...
        'GradObj','off',...
        'GradConstr','off',...
        'Display', 'iter',...
        'Diagnostics', 'off',...
        'FunValCheck', 'off',...
        'TolCon', 1e-3,...
        'TolFun', 1e-5,...
        'TolX', 1e-10,...
        'MaxSQPIter',1000);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Newton configuration
    [X_o_res_vec, E_vec] = combi_optim_02(X_o,phys,s_master,pos, Y_LK_vm, mcmc_options,options, model_fun, constraints_fun);
     
end
save results_snl_w1_11 L_max gdim X_o_res_vec s_master phys pos Y_LK_vm
%% plot the sequence results
plot_time_sequence_pretty(X_o_res_vec(1:end,:),s_master,phys,pos, Y_LK_vm);
    



% % % load X_o_dros_embryo_04.mat;
% % % % % [xc yc zc] = get_xyz_clks(X_o);xc(1) = 0;yc(1) = 0;zc(1) = 0;X_o = 0.36*[xc(:)' yc(:)' zc(:)'];
% % % d = dros_embryo(4,sh_basis(4,120), X_o);plot_pretty(d);view(-90,0);
% % % load ge;
% % % L = 32;
% % % nico = 6;
% % % 
% % % s = sh_surface(L,sh_basis(L,120));
% % % s.xc = sh_surface.tr_xc(gtw,L);s = sh_rot(s,0, 0, pi);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'twist_globe.fig');myaa;print -dtiff -r600 twist_globe.tif;
% % % d = add_sf(d,'twist',s);
% % % 
% % % s.xc = sh_surface.tr_xc(gsnl,L);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'snail_globe.fig');myaa;print -dtiff -r600 snail_globe.tif;
% % % d = add_sf(d,'snail',s);
% % % 
% % % s.xc = sh_surface.tr_xc(ghkb,L);
% % % figure;plot_globe(s,nico);view(-168,10);saveas(gcf,'hkb_globe.fig');myaa;print -dtiff -r600 hkb_globe.tif;
% % % d = add_sf(d,'hkb',s);
% % % 
% % % plot_field(d,nico,1:3);view(24,-60);saveas(gcf,'overlayed.fig');myaa;print -dtiff -r600 overlayed.tif;
% % % plot_field(d,nico,1,1);view(24,-60);saveas(gcf,'twist.fig');myaa;print -dtiff -r600 twist.tif;
% % % plot_field(d,nico,2,2);view(24,-60);saveas(gcf,'snail.fig');myaa;print -dtiff -r600 snail.tif;
% % % plot_field(d,nico,3,3);view(24,-60);saveas(gcf,'hkb.fig');myaa;print -dtiff -r600 hkb.tif;
% % % save d_o d


