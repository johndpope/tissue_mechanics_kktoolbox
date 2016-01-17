function plot_time_sequence_anisotropy(X_o_res_vec,phys)

if nargin <= 5, az = -40;el = 0;end
% close all;
h = dfig;
%set(h,'Color',[0 0 0]);
view(az,el);axis equal;axis off;

%% generate basis for plotting
nicos = 4;      %
[X, C]=surface_mesh.sphere_mesh_gen(nicos);
% [X, C] = sphere_mesh_gen_rand(80);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
L_max = shp_surface.get_L_max(X_o_res_vec(1,:));
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
Y_LK  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;

%% calcualte anisotropy for the initial configuration
X_o = phys.uX_o';
[xclks yclks zclks] = get_xyz_clks(X_o);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
m = surface_mesh(X,C);
m = m.props;

%%% calculate anisotropy for current configuration
X_slave = X_o_res_vec(1,:)';
[xclks yclks zclks] = get_xyz_clks(X_slave);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
mcurr = surface_mesh(X,C);
mcurr = mcurr.props;

anisotropy = (mcurr.quality-m.quality).^2./max((mcurr.quality-m.quality).^2);
%%% plot
clf;
colordef white;
h = figure('Color', [1 1 1]);
set(h,'InvertHardcopy','off');
patch('vertices',X,'faces',C,'FaceColor','flat', 'CData', anisotropy, 'FaceAlpha', 1.0, 'EdgeColor', 'none');  % the middle surface
drawnow
rotate3d on;





