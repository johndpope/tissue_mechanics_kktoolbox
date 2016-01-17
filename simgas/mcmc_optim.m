function [X_o_res_vec, E_vec] = mcmc_optim(X_o,phys,s_master,pos, Y_LK_vm,mcmc_options, model_fun, constraints_fun)
 zbnd = phys.zbnd;   
dfig;
    TolCon = mcmc_options.TolCon;
    maxiter = mcmc_options.maxiter ;
    sigQcon = mcmc_options.sigQcon ;       % when constraints are not satisfied:standard deviation of Gaussian representing probability density for drawing samples
    sigQsat = mcmc_options.sigQsat;       % case when constraints are satisfied
    sigE = mcmc_options.sigE     ;


    fac1 = 1/2/sigE^2;
%     [Et, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(xp,phys);
    X_o_res_vec = [];       % store the resulting(optimal) parameters at each pseudo-time point
    E_vec = [];             % store the shape energy
    for ix = 1:size(pos,1)  % start iterations over pseudo time-points
        disp(['************** position: ' num2str(pos(ix,:))]);
        s_master = position(s_master,pos(ix,:));
        phys.X_o_vm = s_master.X_o;
        [xc yc zc] = get_xyz_clks(phys.X_o_vm);
        phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for master
        phys.cm = sum(phys.X_vm,1)./length(phys.X_vm);
        % start optimization iterations
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
            xp = X_o + sigQ.*randn(numel(X_o),1);
            % check for constraints satisfaction
            fit_parms = parmvec_gen(xp,phys);
            [c ceq] = constraints_fun(fit_parms,phys);   % test the energy calculation
            if any(ceq>TolCon), 
                %disp(['constraint violation: ' num2str(ceq(:)')]);
            else   % i.e. if the constraints are satisfied
                sigQ = sigQsat;
                fit_parms = parmvec_gen(xp,phys);
                [Ep, A, V, v_red,AreaDiff, Eb, EADE, E_stretch, E_shear, T, H] = model_fun(fit_parms,phys);
                P = exp(fac1*(Et-Ep));  % probability ratio used for acceptance test
                if P>= 1 || P>=U(jx),
                    xp = recover_parms(fit_parms,phys);
                    X_o = xp;       %% acceptance, set current vector equal to proposed vector
                    disp(['n: ' num2str(jx) ' t: ' num2str(ix) ' V: ' num2str(V/phys.V_o) '  Acceptance: ' sprintf('%g',P) ' E: ' num2str([Et Ep])]);
                    Et  = Ep;   % set the current energy to be equal to the proposed energy
                    subplot(1,2,1);plot_mcmc(over_Basis(X_o),phys);title(num2str(Et));
                    E_vec = [E_vec Et];
                    subplot(1,2,2);cla;plot(E_vec);drawnow;
                end
            end
        end
        X_slave = over_Basis(X_o);
        X_o_res_vec = [X_o_res_vec;X_slave(:)'];
        
        %%% check for intersection
        [xclks yclks zclks] = get_xyz_clks(X_slave);
        X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
        [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
        [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
        if any(X(:,3)<zbnd), zviolation = 1e0;else zviolation = 0;end
        vec = [res_self res_vm zviolation];
        if any(vec), disp(vec);disp('********** Aborting pseudo-time sequence');break;end
    end