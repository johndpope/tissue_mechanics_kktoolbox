function [C M Z] = kk_cluster(SHP_ri, n)
%%% groups the objects in the different views (length(SHP_ri));
%%% Note: what C tells us is the indices of objects into each SHP_ri view that
%%% correspond to each other, i.e. values in C correspond to bin numbers,
%%% and the position is the index into a cell2mat array (below).
%%% M stores the dissimilarities
%%% Z stores the transformed vectors corresponding to the neighbors lists
%%% Trans stores all transformations corresponding to the bins
if nargin ==1, n = 3;end;



%% build a point-cloud for each object in each view (i.e., prepare geometric local descriptor)
[nlists nobjects viewix objix]= generate_neighbor_lists(SHP_ri, n);

%% Determine the correspondence between every object and all objects in other views using Procrustes analysis
[M Z] = distance_matrix(nlists, nobjects, viewix, objix); 

%%  use hiarchical clustering to group the objects
C = cluster(linkage(squareform(tril(M,-1) + triu(M,1)),'single'),...
    'depth',2, 'cutoff', 1e0, 'criterion', 'inconsistent');
disp(sort(C));
disp(['Number of cluster bins: ' num2str(max(C))]);





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nlists nobjects viewix objix] = generate_neighbor_lists(SHP_ri,n)
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M Z] = distance_matrix(nlists, nobjects, viewix, objix)
M = inf*ones(sum(nobjects), sum(nobjects));    % this is the dissimilarities matrix of among all parametric clouds
Z = cell(sum(nobjects),sum(nobjects));Z{sum(nobjects),sum(nobjects)} = [];
%Trans = cell(sum(nobjects),1);Trans{sum(nobjects)} = [];
gld = [];
for ix = 1:size(M,1),       % loop over the rows of M (which is square anyway)
    v1ix  = viewix(ix);
    Ob1ix = objix(ix);
    n1 = nlists{v1ix};
    X1 = n1{Ob1ix};
    for jx = 1:size(M,2),   % now loop over the columns of M
        v2ix  = viewix(jx);
        Ob2ix = objix(jx);
        if ~(v1ix==v2ix),
            n2 = nlists{v2ix};
            X2 = n2{Ob2ix};
            %% calculate the dissimilarity using procrustes analysis
            [m z t] = procrustes(X1,X2,'reflection',false);
            M(ix,jx) = m;
            Z{ix,jx} = z;   % Z = b*Y*T + c , i.e., z = scale*X2*rotation + translation;
            %Trans{ix,jx} = t;
        end
        %         disp([counter ix jx v1ix Ob1ix v2ix Ob2ix M(ix,jx)]);
    end
end
