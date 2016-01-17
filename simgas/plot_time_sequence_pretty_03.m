function plot_time_sequence_pretty_03(X_o,phys,str)

colordef black;
h = figure('Color', [0 0 0]);
set(h,'InvertHardcopy','off');
axis equal;axis off;

%% generate basis for plotting
nicos = 4;      % number of icosahedron subdivisions used for self-intersection test
[X, C]=surface_mesh.sphere_mesh_gen(nicos);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
L_max = shp_surface.get_L_max(X_o(:));
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
phys.Y_LK_self  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
%% generate color data for the uniform colorbar scaling
LMC = [];   % local mean curvature
X_slave = X_o(:)';
[xclks yclks zclks] = get_xyz_clks(X_slave);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, C, 0);
LMC = H(:);
LMC(:) = LMC(:)/max(LMC(:));
LMC(:) = mat2gray(LMC,[-0.3 0.4]);

%% plot the sequence of shapes
X_slave = X_o(:);
[xclks yclks zclks] = get_xyz_clks(X_slave);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%

subplot(1,3,1);
cla;
%plot_shell(X, C,phys.D);
%     patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 1, 'EdgeColor', 'none');  % the middle surface
% maxH = 0.0115;
% minH = -0.0303;
%H = mat2gray(H,[minH maxH]);
patch('Vertices', X, 'Faces', C,'FaceVertexCData',-H, 'FaceColor', 'interp','EdgeColor', 'none','FaceAlpha',1);
%view(-188,16);
%lighting flat;camlight;
%view(0,0);
%lighting phong;camlight;
%caxis([minH maxH]);
ylim([800 900]);
axis equal;axis on;axis tight;
view(0,0);lighting flat;camlight;drawnow;

subplot(1,3,2);
cla;
%plot_shell(X, C,phys.D);
%     patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 1, 'EdgeColor', 'none');  % the middle surface
maxH = 0.0115;
minH = -0.0303;
%H = mat2gray(H,[minH maxH]);
patch('Vertices', X, 'Faces', C,'FaceVertexCData',-H, 'FaceColor', 'interp','EdgeColor', 'none','FaceAlpha',1);
view(-188,16);
lighting phong;camlight;
view(-90,-90);
lighting phong;camlight;
caxis([minH maxH]);
view(-55,20);
axis equal;axis off;

subplot(1,3,3)
axis off;
text(0,0,str);
drawnow
rotate3d on;