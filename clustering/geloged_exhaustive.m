function [S] = geloged_exhaustive(SHP_ri, n)
%%%% UNFINISHED

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
%% Exhaustive calculation of cost function for all possible configurations 
C = zeros(size(M,1),length(SHP_ri));        % cluster data structure initialization
C(:,1) = 1:size(C,1);                       % each object is it's own cluster at the beginning
orphan_penalty = 1e-20;                         % low values panalize orphans

%% generate vector of legal combinations (without configuration repetitions)
%%% the number of clusters is determined by Stirling number of the second
%%% kind
St = 0;
for m = 3
    N = 4;
    fac = 0;
    for ix = 0:m
        fac = fac+(-1).^(m-ix)*(factorial(m)/factorial(m-ix)/factorial(ix))*(ix)^N;
    end
    St = St + 1/factorial(m)*fac;
end
disp(St)




%% %
function L = cluster_bin_cost(M,L, v)
v = v(v~=0);  % remove zeros
for ix = 1:length(v),
    for jx = ix+1:length(v)
        L = L * exp(-0.5* M(v(ix),v(jx)));
    end
end






























