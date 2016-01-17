function [X_o_res_vec, E_vec] = newton_optim(X_o,phys,s_master,pos, Y_LK_vm,options, model_fun, constraints_fun)

zbnd = phys.zbnd;
%%% start pseudo-time iterations. For each iteration a full optimization is performed
    X_o_res_vec = [];
    E_vec = [];
    exitflag_vec = [];
    
    for ix = 1:size(pos,1)
        s_master = position(s_master,pos(ix,:));
        phys.X_o_vm = s_master.X_o;
        [xc yc zc] = get_xyz_clks(phys.X_o_vm);
        phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for vitteline membrane
        phys.cm = sum(phys.X_vm,1)./length(phys.X_vm);
        %%% optimize configurational energy
        fit_parms = parmvec_gen(X_o,phys);
        [X_result,fval, exitflag, output,lambda_minimization, grad, hessian] = ...
            fmincon(model_fun,fit_parms, [],[],[],[],[],[], constraints_fun,options,phys);
        X_result = recover_parms(X_result, phys);
        X_slave = over_Basis(X_result);
        X_o_res_vec = [X_o_res_vec;X_slave(:)'];
        E_vec = [E_vec fval];
        exitflag_vec = [exitflag_vec exitflag];
        X_o = X_result;     % assign the input to the next iteration as the result of the previous
        
        %% check for intersection
        [xclks yclks zclks] = get_xyz_clks(X_slave);
        X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
        [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
        [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
        if any(X(:,3)<zbnd), zviolation = 1e0;else zviolation = 0;end
        vec = [res_self res_vm zviolation];
        if any(vec), disp(vec);disp('********** Aborting pseudo-time sequence');break;end
    end