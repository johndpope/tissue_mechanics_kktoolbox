function [S,C] = geloged_mcmc_optimization_02(SHP_ri, n)
%%%% Works but slow; Under construction

%%% Performs the generalized local geometric descriptor algorithm on the
%%% SHP shape clouds given in S_cell and uses the optional weighting in Score_cell in the final fusion step.
%%% n is the number of neighbors to consider.
%%% Output: Returns the fused cloud

S = [];
%% build a point-cloud for each object in each view (i.e., prepare geometric local descriptor)
nobjecs = zeros(length(SHP_ri),1);  % vector to store the number of objects for each view
objix   = [];                       % vector to store the index of objects as we go through the views
viewix  = [];                       % vector to store the view index for each object as we go through the views
nlists = {};
for ix = 1:length(SHP_ri)
    T = SHP_ri{ix};
    nobjects(ix) = size(T,1);
    nlist = {}; % neighbor list of each shape
    for row = 1:size(T,1),      % loop over all shapes
        X_o = T(row,:);         % this is our current object
        d = zeros(size(T,1));    % initialize the vector of distances to other objects
        for rix = 1:size(T,1),  % loop over all shapes
            d(rix) = shp_d(T(row), T(rix));        % calculate the distance between shape "row" and shape "rix"
        end
        %%% construct neighbor list
        [d IX] = sort(d);       % sort d in ascending order and obtain the indices of the partners
        d = d(1:n+1);        % take the n first elements of d
        IX = IX(1:n+1);      % indices of the n closes neighbors (plus self)
        neigh = T(IX(2:end),:);    % list of coordinates of the n neighbors
        neigh = neigh-X_o(ones(size(neigh,1),1),:);        % center point cloud around the bead in CLK parameter space
        nlist{row} = neigh;

        viewix = [viewix ix];
        objix  = [objix row];
    end
    nlists{ix} = nlist;
end
%% Determine the correspondence between every object and all objects in other views using Procrustus analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = inf*ones(sum(nobjects), sum(nobjects));    % this is the dissimilarities matrix of among all parametric clouds
gld = [];
counter = 0;
for ix = 1:size(M,1),       % loop over the rows of M (which is square anyway)
    v1ix  = viewix(ix);
    Ob1ix = objix(ix);
    n1 = nlists{v1ix};
    X1 = n1{Ob1ix};
    for jx = 1:size(M,2),   % now loop over the columns of M
        counter = counter + 1;
        v2ix  = viewix(jx);
        Ob2ix = objix(jx);

        if ~(v1ix==v2ix),
            n2 = nlists{v2ix};
            X2 = n2{Ob2ix};
            %% calculate the similarity using procrustes analysis
            M(ix,jx) = procrustes(X1,X2);
        end
        %         disp([counter ix jx v1ix Ob1ix v2ix Ob2ix M(ix,jx)]);
    end
end
%%% note: the matrix M is a distance matrix
S = M;

%% Generate the starting guess and initial quantities for the cluster data structure
C = zeros(size(M,1),length(SHP_ri));        % cluster data structure initialization
C(:,1) = 1:size(C,1);                       % each object is it's own cluster at the beginning
orphan_penalty = 1e-20;                         % low values panalize orphans
v_bins = 5;                                 % guess at the number of bins
accept = 0;
max_iter = 1e6;
accept_vec = zeros(length(max_iter),1);

%% start the mcmc iterations
for jx = 1:max_iter
    %% choose the type of mcmc move (diffusion, birth, death)
% %     mcmc_move = ceil(10.*rand(1,1)); % choose type of move
% %     if mcmc_move < 7, mcmc_move = 1;
% %     elseif mcmc_move<9,mcmc_move = 2;
% %     else mcmc_move = 3;end
mcmc_move = 1;
    switch(mcmc_move)
        case 1          %% Diffusion
            %% pick random entry in C switch with another random entry and evaluate
            %% acceptance probability
            cargo_r = ceil(size(C,1).*rand(1,1)); % choose cluster bin index
            c = C(cargo_r,:);
            if sum(c(:)) >0     % then we have an occupied (legal) cluster bin
                cargo_c = ceil(length(SHP_ri).*rand(1,1)); % choose which entry to change
                cargo_ix = C(cargo_r, cargo_c);     % this is the index to be removed from here and placed somewhere else
                target_r = ceil(size(C,1).*rand(1,1)); % choose cluster bin index
                t = C(target_r,:);
                if ~(target_r==cargo_r)
                    %if sum(t(:)) >0
                    target_c = ceil(length(SHP_ri).*rand(1,1)); % choose which entry to change
                    target_ix = C(target_r, target_c);     % this is the index to be replaced and needs to go somethere else
                    if ~(target_ix ==0 && cargo_ix ==0)
                        %disp([cargo_r cargo_c cargo_ix target_r target_c target_ix]);
                        %%%
                        L_before = 1;
                        if sum(c==0)==length(SHP_ri)-1, L_before = L_before * orphan_penalty;   % i.e. if it is an orphan in a bin
                        elseif sum(c==0) == length(SHP_ri),     % do nothing
                        else
                            L_before = cluster_bin_cost(M,L_before, c);
                        end

                        %if L_before==0, disp([cargo_r cargo_c cargo_ix target_r target_c target_ix]); warning('L_before is zero');end

                        if sum(t==0)==length(SHP_ri)-1, L_before = L_before * orphan_penalty; % i.e. if it is an orphan
                        elseif sum(t==0) == length(SHP_ri),     % do nothing
                        else L_before = cluster_bin_cost(M,L_before, t); end

                        %if L_before==0, disp([cargo_r cargo_c cargo_ix target_r target_c target_ix]); warning('L_before is zero');end
                        %%%
                        c = C(cargo_r,:);
                        t = C(target_r,:);
                        c(cargo_c) = t(target_c);
                        t(target_c) = cargo_ix;
                        %%%
                        L_after = 1;
                        if sum(c==0)==length(SHP_ri)-1, L_after = L_after * orphan_penalty;   % i.e. if it is an orphan
                        elseif sum(c==0) == length(SHP_ri),     % do nothing
                        else
                            L_after = cluster_bin_cost(M,L_after, c);
                        end
                        %if L_after==0, warning('L_after is zero');end
                        if sum(t==0)==length(SHP_ri)-1, L_after = L_after * orphan_penalty; % i.e. if it is an orphan
                        elseif sum(t==0) == length(SHP_ri),     % do nothing
                        else
                            L_after = cluster_bin_cost(M,L_after, t);
                        end
                        %if L_after==0, warning('L_after is zero');end

                        if L_before == 0, alpha = 1;
                        else alpha = min([1 L_after/L_before]);
                        end

                        if alpha==1, accept=1;
                        elseif alpha<rand(1), accept = 0;
                        end

                        if accept, % update the Cluster data-structure
                            temp = cargo_ix;
                            C(cargo_r,cargo_c) = target_ix;
                            C(target_r,target_c) = temp;
                        end
                        %disp([alpha L_after L_before accept]);
                        %                 disp([alpha L_after/L_before accept]);
                        %                 disp(C);
                        %                 pause(0.2);
                        accept_vec(jx) = accept;
                    end
                    %else accept_vec(jx) = 0;
                    %end
                end
            else accept_vec(jx) = 0;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 2          %% Birth
            cargo_r = ceil(size(C,1).*rand(1,1)); % choose cluster bin index
            c = C(cargo_r,:);
            if sum(c(:)) >0     % then we have an occupied (legal) cluster bin
                cargo_c = ceil(length(SHP_ri).*rand(1,1)); % choose which entry to remove
                cargo_ix = C(cargo_r, cargo_c);     % this is the index to be removed from here and placed somewhere else
            end
            L_before = 1;
            if sum(c==0)==length(SHP_ri)-1, L_before = L_before * orphan_penalty;   % i.e. if it is an orphan in a bin
            elseif sum(c==0) == length(SHP_ri),     % do nothing, it is an empty bin
            else
                L_before = cluster_bin_cost(M,L_before, c);
            end
            
            L_after = 1;
            if sum(c==0)==length(SHP_ri)-1, L_after = L_after * orphan_penalty;   % i.e. if it is an orphan
            elseif sum(c==0) == length(SHP_ri),     % do nothing
            else
                L_after = cluster_bin_cost(M,L_after, c);
            end
            
            
            if L_before == 0, alpha = 1;
            else alpha = min([1 L_after/L_before* v_bins/(size(C,1) + 1)* orphan_penalty]);
            end
            if alpha==1, accept=1;
            elseif alpha<rand(1), accept = 0;
            end
            if accept, % update the Cluster data-structure
                t = zeros(length(SHP_ri),1);
                t(1) = cargo_ix;
                C(cargo_r,cargo_c) = 0;
                C(end+1,:) = t;     % add the additional bin
            end
            accept_vec(jx) = accept;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3          %% Death
            indx = find(sum(C,2)==0); %#ok<EFIND>
            if ~isempty(indx), 
                alpha = min([1 v_bins/(size(C,1) - 1)]);
                if alpha==1, accept=1;elseif alpha<rand(1), accept = 0; end
                if accept, C(indx(1),:) = [];end
                accept_vec(jx) = accept;
            end
    end

end
disp(C);
disp('Occupied bins'); 
disp(length(C));
C(find(sum(C,2)==0),:) = [];
disp(C);

itercount = 1:length(accept_vec);
figure;plot(cumsum(accept_vec)./itercount * 100);







%% %
function L = cluster_bin_cost(M,L, v)
v = v(v~=0);  % remove zeros
for ix = 1:length(v),
    for jx = ix+1:length(v)
        L = L * exp(-0.5* M(v(ix),v(jx)));
    end
end






























