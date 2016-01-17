function plot_shell(X,C, d)
% the usual X and C, and the thickness d
%% plot the slave surface
cla;

% patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 0.8);  % the middle surface


%%% now that we have the normals at each point we calculate the surface that 
%%% is a distance d removed from it
[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, C, 0);
X_in = X + N .* d/2;  %%% now generate the inner surface
patch('vertices',X_in,'faces',C,'FaceColor','blue', 'FaceAlpha', 0.8);
X_out = X - N .* d/2;  %%% now generate the outer surface
patch('vertices',X_out,'faces',C,'FaceColor','blue', 'FaceAlpha', 0.8);





