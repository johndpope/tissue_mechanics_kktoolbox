function X_o = phase_diagram(OM,A_o_vec, q,v_o_vec, a_o_bar_vec, X_o,...
    Y_LK_tri,t, p,  flag,lambda_shear, filename, mc_opts, options)
% global r_o v_o da_o gdim p t L K wg q
global itercount
counter = 1;
gdim = size(Y_LK_tri,1);

for step = 1:length(a_o_bar_vec)
    v_o = v_o_vec(step); a_o_bar = a_o_bar_vec(step); itercount = 0;
    disp(['Starting minimization number ' num2str(counter)]);
    if length(A_o_vec)>1,A_o = A_o_vec(step);else A_o = A_o_vec;end
    v_red_o = v_o./(4/3*pi*(A_o/4/pi).^(3/2));
    disp([A_o v_o v_red_o a_o_bar]);

    % %%%%%%%%%%%%%%%%% Start minimization SQP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OM == 1,

        [X_o,fval, exitflag, output,lambda_minimization, grad, hessian] = ...
            fmincon(@objfun_ade,X_o,[],[],[],[],[],[], @confuneq_ade,options,A_o,...
            q,Y_LK_tri,v_o,a_o_bar, lambda_shear,flag);
        [E v da] = objfun_ade(X_o,A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear, 1);
%         confuneq_ade(X_o,A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear, 1);
        str = sprintf('Bending energy = %.6f   v = %.4f v_o = %.4f  da = %.4f da_o = %.4f',full(fval),full(v), full(v_o), full(da), full(a_o_bar_vec(step)));disp(str);
        disp(lambda_minimization.eqnonlin);
        lambda_vec(counter) = lambda_minimization;
        save last_X_o X_o
        nc = round(length(X_o)/3);
        xclks = X_o(1:nc);yclks = X_o(nc+1:2*nc); zclks = X_o(2*nc+1:3*nc);
        lmax = sqrt(length(X_o)/3)-1;
        str = sprintf('save X_o_L%d_A_%.1f_V%.1f_da%.5f_lambda_shear_%.2f_%s.mat X_o fval exitflag output lambda_minimization grad hessian;',lmax, A_o,v_o,a_o_bar, lambda_shear,filename);eval(str);
        X_o_trs = trs_invariance_gen(X_o);
        str = sprintf('save X_o_trs_L%d_A_%.1f_V%.1f_da%.5f_lambda_shear_%.2f_%s.mat X_o_trs;',lmax, A_o,v_o,a_o_bar, lambda_shear,filename);eval(str);
        v_o_rec(counter) = v_o; a_o_bar_rec(counter) = a_o_bar;wb_rec(counter) = E; q_rec(counter) = q;
        x_vec(counter,:) = X_o;%fix(clock)
        counter = counter + 1;%str = sprintf('save %s;',filename);eval(str);

    elseif OM ==2         %%% Start minimization MC
        clks = X_o;
        X_o = pd_montecarlo(mc_opts,...  % 1 = precalcualte constraints-satisfying shapes
            A_o, 2/pi, v_o, a_o_bar,X_o, Y_LK_tri,t, p,  0, lambda_shear,'_oblateMS.mat');
        if isempty(X_o), X_o = clks; disp('No shape satisfied constrains. Returning same shape.');return;end;
        clks = X_o;
        [E v A] = objfun_ade(X_o,A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear, 2);

        str = sprintf('ADE energy = %.6f   RedVol = %.4f',full(E),full(v));disp(str);
        save last_X_o X_o ;
        lmax = sqrt(length(X_o)/3)-1;
        str = sprintf('save X_o_MC_L%d_A_%.1f_V%.1f_da%.5f_%s.mat X_o;',lmax, A_o,v_o,a_o_bar, filename);eval(str);
%         X_o_trs = trs_invariance_gen(X_o);
%         str = sprintf('save X_o_trs_MC_L%d_A_%.1f_V%.1f_da%.5f_%s.mat X_o_trs;',lmax, A_o,v_o,a_o_bar, filename);eval(str);
    end
end
