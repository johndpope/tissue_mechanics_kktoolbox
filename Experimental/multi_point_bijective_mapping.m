load XF_cube
X = X';
F = F';
T = triangulation(F,X);

ixN = 1;
ixS = size(X,1);

%[t,p,~,~] = surface_mesh.bijective_map_gen(X, F, L, 1); % %%%% Generate Bijective Map
m = surface_mesh(X,F);
m = m.edge_info;
edg = unique(m.E,'rows');
A = sparse(size(X,1));
for eix = 1:size(m.E,1)
    A(m.E(eix,1), m.E(eix,2)) = 1;
end
A(1,1) = 0;
%% construct the graph
G = graph(A);
% find links
L = {};
for ix = 1:size(X,1)
    nbrs = neighbors(G,ix);
    L{ix} = nbrs';
end

[t, A, b] = surface_mesh.latitude_calc(L, ixN, ixS); %save data_latitude    % Calculate theta (latitude values associated with each vertex)
figure;patch('Vertices',X,'Faces',F,'FaceVertexCData',t, ...
              'FaceColor','interp', 'EdgeColor','k');
axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;

[p, A, b, dtline, W] = surface_mesh.longitude_calc(X(:,1), X(:,2), X(:,3),...
                                                   t, A, F, L, ixN, ixS); %save data_longitude
