function P = minimize_translation_global(Po,Xcell,parmix)
Xcell_try = Xcell;
%res = objective_LM(Po);

[P,resnorm,residual] = lsqnonlin(@objective_LM,Po,[],[],...
    optimset('Algorithm', 'levenberg-marquardt',...
    'LargeScale','off', ...
    'MaxFunEvals', 1000,...
    'Display', 'off'));

    function res = objective_LM(Pall)
        verbose = 0;
        res = [];
        
        for ix = 1:size(parmix,1),   % loop over the 3 sets of parameters
            P = Pall((3*(ix-1))+1:3*ix);
            set1ix = parmix(ix,1);
            Xcell_try{set1ix} = translate_rotate_XF(Xcell{set1ix}, [P(:)' 0 0 0]); % only do translation
            set2ix = parmix(ix,2);
            Xcell_try{set2ix} = translate_rotate_XF(Xcell{set2ix}, [P(:)' 0 0 0]); % only do translation
        end
        if verbose, figure(1);clf;end
        for ix = 1:size(Xcell,1),       % loop over the data pairs
            %%% calculate residual vector given the current Xcell_try
            res  = [res dist_mx_calc(Xcell_try{ix,1},Xcell_try{ix,2})];%% calculate the distance matrix
            if verbose
                Xtr = Xcell_try{ix,1};
                X   = Xcell_try{ix,2};
                plot3(X(:,1),X(:,2),X(:,3),'*b');hold('on');
                plot3(Xtr(:,1),Xtr(:,2),Xtr(:,3),'*r');hold on;
                xlabel('x');ylabel('y');zlabel('z');view(2);axis equal;drawnow;
            end
        end
        if verbose, hold off;end
    end

% % % % % % %%% now refine with simplex
% % disp('Starting simplex refinement');
% % [P,fval] = fminsearch(@objective_simplex, P, ...
% %     optimset('Display','iter', 'MaxIter', 2000, 'TolX', 1e-3 ));
% % 
% %     function res = objective_simplex(Pall)
% %         verbose = 1;
% %         res = [];
% %         
% %         for ix = 1:size(parmix,1),   % loop over the 3 sets of parameters
% %             P = Pall((3*(ix-1))+1:3*ix);
% %             set1ix = parmix(ix,1);
% %             Xcell_try{set1ix} = translate_rotate_XF(Xcell{set1ix}, [P(:)' 0 0 0]); % only do translation
% %             set2ix = parmix(ix,2);
% %             Xcell_try{set2ix} = translate_rotate_XF(Xcell{set2ix}, [P(:)' 0 0 0]); % only do translation
% %         end
% %         if verbose, figure(1);clf;end
% %         for ix = 1:size(Xcell,1),       % loop over the data pairs
% %             %%% calculate residual vector given the current Xcell_try
% %             res  = [res dist_mx_calc(Xcell_try{ix,1},Xcell_try{ix,2})];%% calculate the distance matrix
% %             if verbose
% %                 Xtr = Xcell_try{ix,1};
% %                 X   = Xcell_try{ix,2};
% %                 plot3(X(:,1),X(:,2),X(:,3),'*b');hold('on');
% %                 plot3(Xtr(:,1),Xtr(:,2),Xtr(:,3),'*r');hold on;
% %                 xlabel('x');ylabel('y');zlabel('z');view(2);axis equal;drawnow;
% %             end
% %         end
% %         if verbose, hold off;end
% %         res = sum(res.^2);
% % 
% %     end
end

function res = dist_mx_calc(X,Xtr)
C0x = X(:,1);C0x = C0x(:,ones(size(Xtr,1),1));
C0y = X(:,2);C0y = C0y(:,ones(size(Xtr,1),1));
C0z = X(:,3);C0z = C0z(:,ones(size(Xtr,1),1));

Cdrx = Xtr(:,1)';Cdrx = Cdrx(ones(size(X,1),1), :);
Cdry = Xtr(:,2)';Cdry = Cdry(ones(size(X,1),1), :);
Cdrz = Xtr(:,3)';Cdrz = Cdrz(ones(size(X,1),1), :);

Rsq = ((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances squared
res = min(Rsq);  % row vector of  minimum distances squared
end