function [S] = geloged_deterministic_annealing(SHP_ri, n)
%%% Generates for each object, a vector comprised of its own coefficients
%%% and n closest neighbor coefficients. In this high dimensional space the
%%% deterministic annealing takes place for clustering to objects
%%% Output: Returns the fused cloud
S = [];
%% build a point-cloud for each object in each view (i.e., prepare local morphological context)
nobjecs = zeros(length(SHP_ri),1);  % vector to store the number of objects for each view
objix   = [];                       % vector to store the index of objects as we go through the views
viewix  = [];                       % vector to store the view index for each object as we go through the views
V       = [];                       % matrix to with each column representing one feature vector
for ix = 1:length(SHP_ri)       % loop over all views
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
        [d IX] = sort(d);    % sort d in ascending order and obtain the indices of the partners
        IX = IX(1:n+1);      % indices of the closest neighbors (plus self)
        neigh = T(IX(2:end),:);    % coefficients (coordinates) of the n neighbors
        v = [X_o(:);neigh(:)];
        V = [V v];
        viewix = [viewix ix];
        objix  = [objix row];
    end
    nlists{ix} = nlist;
end
%% %%%%%%%%%%%%%%%%%%%%%%%
disp(V);
S = V;
%% Deterministic annealing is now required to generate the clusters

y = sum(V,2)./size(V,2); % Starting value for the average coefficient vector
















