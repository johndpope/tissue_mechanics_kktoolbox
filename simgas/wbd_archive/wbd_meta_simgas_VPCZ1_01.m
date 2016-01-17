%%%% modeling ventral furrow invagination
clc;close all;global counter;counter = 0;
in.runName = sprintf('DE_08_VPCZ1_00_%2f',rand);
%% Step 0: configure L_max  gdim and optimization
in.L_max = 6;           % maximum order of the spherical harmonics parameterization
in.L_max_gef = 5;
in.gdim = 30;
in.nicos = 2;       % for self intersection test (number of icosahedron subdivisions)
in.nicovm = 2;
in.optim_method = 2;
in.pos = [ 0 0 0];
in.model_fun = @DG_energy_shell_08_production;
in.constraints_fun = @DG_constraints;
%% step 1: define data file names for morphologies and precalculated gene expression and intersection information
in.fn_master = 'X_o_dros_embryo_04.mat';
in.fn_slave  = 'X_o_dros_embryo_04.mat';
in.fn_start  = 'X_o_dros_embryo_04_ind1.mat';
in.fn_ge     = 'ge.mat';
in.fn_data_self = 'DATA_self_nico_2.mat';
in.fn_tdb_data  = 'TDB_nico2.mat';
%% configure geometrical and physical quantities -- check units!! %%%%%
phys.D  = 15;  % units micron
DEL_el_sq = 0.5^2;%0.28^2; %large values make bending harder, low make shear harder
% low values make shear and stretch harder
phys.kb = 1; %units
phys.k_stretch = phys.kb/DEL_el_sq;
phys.k_shear = phys.k_stretch/2;
phys.a3 = -2;
phys.a4 = 8;
phys.b1 = 0.7;
phys.b2 = 0.75;

% phys.Young = 1000;
% phys.Poiss = 0.3;
% phys.miu = phys.Young/2/(1+phys.Poiss);
% phys.lambda = phys.Poiss*phys.Young/(1-2*phys.Poiss)/(1+phys.Poiss);

phys.lambda_ge = 2.0e0;  % scale for the gene expression field
in.phys = phys;
%% %% configure MCMC method
in.mcmc_options.TolCon          = 1e-3;
in.mcmc_options.maxiter         = 5e3;
in.mcmc_options.sigQcon         = 1e-2;       % when constraints are not satisfied:standard deviation of Gaussian representing probability density for drawing samples
in.mcmc_options.sigQsat         = 1e-2;       % case when constraints are satisfied
in.mcmc_options.sigE            = 5.0e-4;       % standard deviation of the energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of MCMC configuration
%% %% configure matlab's SQP
in.newton_options =   optimset(...%'plotfcn',@DG_opt_plot,...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% configure a vector of 'in' structs and submit in a loop if desired
DEL_el_sq_vec = [1.0:-0.2:0.1].^2;       % provide a range for DEL_el_sq
in_vec = cell(numel(DEL_el_sq_vec),1);
runName = in.runName;
for ix = 1:numel(DEL_el_sq_vec)
    disp(['------------------' num2str(ix) ' of ' num2str(numel(DEL_el_sq_vec)) ' ----------']);
    in.DEL_el_sq = DEL_el_sq_vec(ix);
    in.runName = [runName 'DEL_' num2str(sqrt(in.DEL_el_sq_vec(ix)))];
    [X_o_res_vec, E_vec, exitflag_vec] = submit_DE(in);%% submit
end

%% optional plotting
plot_time_sequence_pretty_02(X_o_res_vec(1:end,:),s_master,phys,pos, Y_LK_vm);
plot_time_sequence_pretty(X_o_res_vec(1:end,:),s_master,phys,pos, Y_LK_vm);