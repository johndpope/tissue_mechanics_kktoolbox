function [S] = geloged(SHP_ri, n)
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
M = ones(sum(nobjects), sum(nobjects));    % this is the dissimilarities matrix of among all parametric clouds
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
        disp([counter ix jx v1ix Ob1ix v2ix Ob2ix M(ix,jx)]);
    end
end
%%% note: the matrix M is a distance matrix
S = M;

%% Obtain a starting guess for the adjacency matrix
A = zeros(size(M));
for ix = 1:size(M,1),
    [B indx] = sort(M(ix,:));
    A(ix,indx(1:n-1)) = 1;
end
disp(A);
cost = sum(M(A==1)); disp(cost);

edge_list = generate_edge_list(A);disp(edge_list);
%% 
% % % % edge_list = generate_edge_list(A);disp(edge_list);
% % % % usnet = edge_list2net(edge_list); % format the edgelist for the loop counting process
% % % % net = sort_net(usnet);
% % % % num_nodes = length(net);
% % % % num_edges = calc_num_edges(net);
% % % % disp(['  Net:  nodes = ' num2str(num_nodes) ' edges = ' num2str(num_edges)]);
% % % % figure; plot_net(net); % plot of the network
% % % % title('Network');
% % % % 
% % % % % % net = reduce_net(net);     % removes all 1-connected nodes and their corresponding edges
% % % % % % num_nodes = length(net);
% % % % % % num_edges = calc_num_edges(net);
% % % % % % disp(['   Reduced net:  nodes = ' num2str(num_nodes) ' edges = ' num2str(num_edges)]);
% % % % 
% % % % % stn = get_starting_node(net);         % give the path a nearly optimal starting node
% % % % stn = 1;
% % % % path = net(stn).node;               % initialize the path
% % % % current_edge = net(stn).edges(1);   % initialize the first edge
% % % % loop_list = [];                   % initialize the loop list
% % % % iterations = 0;              % initialize the number of algorithm steps
% % % % 
% % % % while (length(path)>1 || ~isempty(current_edge))
% % % %     [net,path,current_edge,loop_list] = iterate_tree(net,path,current_edge,loop_list);
% % % %     iterations = iterations+1;
% % % % end
% % % % num_loops = length(loop_list);
% % % % disp(['    It took ' num2str(iterations) ' steps to complete the ILCA']);
% % % % disp(['     There are ' num2str(num_loops) ' loops in the net']);


%%
function E = generate_edge_list(A)
[r c] = ind2sub(size(A), find(A==1));
E = [r(:) c(:)];






























