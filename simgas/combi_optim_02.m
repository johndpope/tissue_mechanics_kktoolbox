function [X_o_res_vec, E_vec, exitflag_vec] = combi_optim_02(X_o,phys,s_master,pos, Y_LK_vm,mcmc_options, options, model_fun, constraints_fun)
plot_flag = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%% First perform MCMC to obtain configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%% that satisfy the constraints
dfig;
TolCon = mcmc_options.TolCon;
maxiter = mcmc_options.maxiter ;
sigQcon = mcmc_options.sigQcon ;       % when constraints are not satisfied:standard deviation of Gaussian representing probability density for drawing samples
sigQsat = mcmc_options.sigQsat;       % case when constraints are satisfied
sigE = mcmc_options.sigE     ;
sigQ = sigQsat;

fac1 = 1/2/sigE^2;
%     [Et, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(xp,phys);
X_o_res_vec = [];       % store the resulting(optimal) parameters at each pseudo-time point
E_vec = [];             % store the shape energy
E_vec_mcmc = [];
exitflag_vec = [];


% start iterations over pseudo time-points
for ix = 1:size(pos,1)  
%     disp(['************** position: ' num2str(pos(ix,:))]);
%     s_master = position(s_master,pos(ix,:));
%     phys.X_o_vm = s_master.X_o;
%     [xc yc zc] = get_xyz_clks(phys.X_o_vm);
%     phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for master
%     phys.cm = sum(phys.X_vm,1)./length(phys.X_vm);
    
    %%    start optimization iterations
    fit_parms = parmvec_gen(X_o,phys);
    [c ceq] = constraints_fun(fit_parms,phys);   % test the constraints calculation
    if any(ceq)>TolCon
        Et = inf;
        sigQ = sigQcon;
        disp('Starting config. does not satisfy constraints... using Et = inf');
    else
        [Et, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);
    end
    U = rand(maxiter,1);    % random number (0-1) drawn from uniform distribution
    for jx = 1:maxiter      % run the optimization for a particular time point X_o
        % generate a new sample
        xp = X_o(:) + sigQ.*randn(numel(X_o),1);
        % check for constraints satisfaction
        fit_parms = parmvec_gen(xp,phys);
        [c ceq] = constraints_fun(fit_parms,phys);   % test the energy calculation
        if any(ceq>TolCon),
            %if plot_flag, subplot(1,2,1);plot_mcmc(over_Basis(recover_parms(fit_parms,phys)),phys);title(num2str(Et));end
            %disp(['constraint violation: ' num2str(ceq(:)')]);
        else   % i.e. if the constraints are satisfied
            sigQ = sigQsat;
            fit_parms = parmvec_gen(xp,phys);
            [Ep, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);
            %disp([Ep Eb E_stretch E_shear]);
            P = exp(fac1*(Et-Ep));  % probability ratio used for acceptance test
            if P>= 1 || P>=U(jx),
                xp = recover_parms(fit_parms,phys);
                X_o = xp;       %% acceptance, set current vector equal to proposed vector
                disp(['n: ' num2str(jx) ' t: ' num2str(ix) ' V: ' num2str(V/phys.V_o) '  Acceptance: ' sprintf('%g',P) ' E: ' num2str([Et Ep])]);
                Et  = Ep;   % set the current energy to be equal to the proposed energy
                if plot_flag, subplot(1,2,1);plot_mcmc(over_Basis(X_o),phys);title(num2str(Et));end
                E_vec_mcmc = [E_vec_mcmc Et];
                if plot_flag,
                    subplot(1,2,2);cla;plot(E_vec_mcmc);
                    str = sprintf('Eb: %.2f  Esh: %.2f   Est: %.2f', Eb/Eb,full(E_shear)/Eb, full(E_stretch)/Eb);
                    title(str);
                    drawnow;
                end
            end
        end
    end
    
    X_slave = over_Basis(X_o);
    %X_o_res_vec = [X_o_res_vec;X_slave(:)'];
    
    %%% check for intersection
    [xclks yclks zclks] = get_xyz_clks(X_slave);
    X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
    [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
    [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
    vec = [res_self res_vm];
    if any(vec), disp(vec);disp('********** Aborting pseudo-time sequence');break;end
    
    %% %% refine the surface with a Newton method
    fit_parms = parmvec_gen(times_basis(X_slave),phys);    % assign input config. as output from mcmc for that time point
    [X_result,fval, exitflag, output,lambda_minimization, grad, hessian] = ...
        fmincon(model_fun,fit_parms, [],[],[],[],[],[], constraints_fun,options,phys);
    X_result = recover_parms(X_result, phys);
    X_slave = over_Basis(X_result);
    X_o_res_vec = [X_o_res_vec;X_slave(:)'];
    E_vec = [E_vec fval];
    exitflag_vec = [exitflag_vec exitflag];
    %X_o = X_result;     % assign the input to the next iteration as the result of the previous
    %%%%
    
end

% % % % %% %%%%%%%%%%%%%%%%%%%%%%%%% Second perform Newton to obtain refined
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%% configurations
% % % % X_o_res_vec2 = [];
% % % % E_vec2 = [];
% % % % exitflag_vec = [];
% % % %
% % % % for ix = 1:size(X_o_res_vec,1)
% % % %     s_master = position(s_master,pos(ix,:));
% % % %     phys.X_o_vm = s_master.X_o;
% % % %     [xc yc zc] = get_xyz_clks(phys.X_o_vm);
% % % %     phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for vitteline membrane
% % % %     phys.cm = sum(phys.X_vm,1)./length(phys.X_vm);
% % % %     %%% optimize configurational energy
% % % %     fit_parms = parmvec_gen(times_Basis(X_o_res_vec(ix,:)),phys);    % assign input config. as output from mcmc for that time point
% % % %     [X_result,fval, exitflag, output,lambda_minimization, grad, hessian] = ...
% % % %         fmincon(model_fun,fit_parms, [],[],[],[],[],[], constraints_fun,options,phys);
% % % %     X_result = recover_parms(X_result, phys);
% % % %     X_slave = over_Basis(X_result);
% % % %     X_o_res_vec2 = [X_o_res_vec;X_slave(:)'];
% % % %     E_vec2 = [E_vec fval];
% % % %     exitflag_vec = [exitflag_vec exitflag];
% % % %     %X_o = X_result;     % assign the input to the next iteration as the result of the previous
% % % %
% % % %     %% check for intersection
% % % %     [xclks yclks zclks] = get_xyz_clks(X_slave);
% % % %     X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
% % % %     [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
% % % %     [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
% % % %     if any(X(:,3)<zbnd), zviolation = 1e0;else zviolation = 0;end
% % % %     vec = [res_self res_vm zviolation];
% % % %     if any(vec), disp(vec);disp('********** Aborting pseudo-time sequence');break;end
% % % % end
% % % %
% % % %
