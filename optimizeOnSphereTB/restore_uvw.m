function [X_new, norm_C] = restore_uvw(constrainsfun, X, miu, tol,iter,flag, plotflag, az, el)
%%% 
global nvert F
[c, C, gc, A] = constrainsfun(X,0); % we will need the constraints as well as the matrix A (Jacobian of the constraints)
if any(c>=0)
    disp(['Inequality constraint activated: c positive iteration ' num2str(iter)]);act_con = gc(:,c>=0); A = [A act_con];C = [C; c(c>=0)];
elseif any(abs(imag(c))>tol),
    disp(['Inequality constraint activated: c imaginary iteration ' num2str(iter)]);act_con = gc(:,abs(imag(c))>=tol);
    A = [A act_con];C = [C; c(abs(imag(c))>=tol)];
end
A = sparse(A);%if any(abs(imag(reshape(A,size(A,1)*size(A,2),1)))>tol),error('A has imaginary elements');end
norm_C = norm(C);
[L1,U1] = luinc(A'*A,1e-11); [dv, flag1, relres1, iter1, resvec1] = pcg(A'*A,-C, 1e-6, 50, L1, U1);


if flag ==1,       % do a line search
    [c ceq] = constrainsfun(X,0);    nc_old = sum((c>=0));       % the number of active inequality constraints
    for lsix = 1:100, miu = miu-miu/2;   [c ceq] = constrainsfun(X+miu*(A*dv),0); nc_new = sum(any(c>=0));nceq_new = sum(any(ceq>=0));
        if  nc_new <= nc_old,  break;end;
    end 
end
X_new = X + miu*A*(dv);
if any(abs(imag(X_new))>tol),error('X is imaginary');end


%%%% plotting
if plotflag==1,
    u = X(1:(nvert)); v = X(nvert+1:2*nvert); w = X(2*nvert+1:end);
    plot_state(u,v,w,F, iter,0,  az, el);%%% Let's look at the current state
end
if plotflag ==2     % shows the constraints 
    u = X(1:(nvert)); v = X(nvert+1:2*nvert); w = X(2*nvert+1:end);
    subplot(1,2,1);plot_state(u,v,w,F, iter,0,  az, el);%%% Let's look at the current state
    subplot(1,2,2);plot(C,'r.');hold on;
    %axis([0 .4 0 .4]);
    ylim([0 6*pi/length(F)]);
    plot(0 * ones(size(C)),'b-');
    
    hold off
end
