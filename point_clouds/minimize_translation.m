function P = minimize_translation(Po,X,Xdata)

res = objective_LM(Po);

[P,resnorm,residual] = lsqnonlin(@objective_LM,Po,[],[],...
    optimset('Algorithm', 'levenberg-marquardt',...
    'LargeScale','off', ...
    'MaxFunEvals', 1000,...
    'Display', 'on'));

    function res = objective_LM(P)
        
        % objective function for registration minimization
        verbose = 1;
        Xtr = translate_rotate_XF(Xdata, [P(:)' 0 0 0]); % only do translation
        
        %% calculate the distance matrix
        C0x = X(:,1);C0x = C0x(:,ones(size(Xtr,1),1));
        C0y = X(:,2);C0y = C0y(:,ones(size(Xtr,1),1));
        C0z = X(:,3);C0z = C0z(:,ones(size(Xtr,1),1));
        
        Cdrx = Xtr(:,1)';Cdrx = Cdrx(ones(size(X,1),1), :);
        Cdry = Xtr(:,2)';Cdry = Cdry(ones(size(X,1),1), :);
        Cdrz = Xtr(:,3)';Cdrz = Cdrz(ones(size(X,1),1), :);
        
        Rsq = ((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances squared
        md = min(Rsq);  % row vector of  minimum distances squared
        res  = md;
        
        if verbose
            figure(1);clf;
            plot3(X(:,1),X(:,2),X(:,3),'*b');hold('on');
            plot3(Xtr(:,1),Xtr(:,2),Xtr(:,3),'*r');hold off;
            xlabel('x');ylabel('y');zlabel('z');view(2);axis equal;drawnow;
        end
    end

% %%% now refine with simplex
[P,fval] = fminsearch(@objective_simplex, P, ...
    optimset('Display','iter', 'MaxIter', 2000, 'TolX', 1e-3 ));

    function res = objective_simplex(P)
        % objective function for registration minimization
        verbose = 1;
        Xtr = translate_rotate_XF(Xdata, [P(:)' 0 0 0]);
        
        %% calculate the distance matrix
        C0x = X(:,1);C0x = C0x(:,ones(size(Xtr,1),1));
        C0y = X(:,2);C0y = C0y(:,ones(size(Xtr,1),1));
        C0z = X(:,3);C0z = C0z(:,ones(size(Xtr,1),1));
        
        Cdrx = Xtr(:,1)';Cdrx = Cdrx(ones(size(X,1),1), :);
        Cdry = Xtr(:,2)';Cdry = Cdry(ones(size(X,1),1), :);
        Cdrz = Xtr(:,3)';Cdrz = Cdrz(ones(size(X,1),1), :);
        
        Rsq = ((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances squared
        md = min(Rsq);  % row vector of  minimum distances squared
        res  = sum(md(:));
        
        if verbose
            figure(1);clf
            plot3(X(:,1),X(:,2),X(:,3),'*b');hold on
            plot3(Xtr(:,1),Xtr(:,2),Xtr(:,3),'*r');hold off
            xlabel('x');ylabel('y');zlabel('z');view(2);axis equal;drawnow
        end
        
    end
end