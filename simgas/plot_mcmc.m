function plot_mcmc(X_o, phys)

cla;
ys = 700;
ye = 900;
% % shp_slice_x(X_o, 100);

%% plot the slave surface
[xclks yclks zclks] = get_xyz_clks(X_o);
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
% view(-46,8);axis equal;axis on;
% xlim([-400 400]);ylim([-400 400]);zlim([-200 200]);
X(X(:,2)<ys,:) = [];
X(X(:,2)>ye,:) = [];
kk_plot3(X);

%patch('vertices',X,'faces',phys.C_slave,'FaceColor','red', 'FaceAlpha', 0.8);


% uncomment below to also display the thickness surface
% % %%% now that we have the normals at each point we calculate the surface that 
% % %%% is a distance d removed from it
% % [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, phys.C_slave, 0);
% % X_in = X + N .* phys.D;  %%% now generate the inner surface
% % patch('vertices',X_in,'faces',phys.C_slave,'FaceColor','blue', 'FaceAlpha', 0.8);
% % hold on


 % ucomment below to also plot the vitteline (master) surface
% % %% plot master
% % patch('vertices',phys.X_vm,'faces',phys.C_vm,'FaceColor','green', 'FaceAlpha', 1);

axis equal;axis off;axis tight;view(0,0);lighting flat;camlight;drawnow;