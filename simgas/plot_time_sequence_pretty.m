function plot_time_sequence_pretty(X_o_res_vec,s_master,phys,pos, Y_LK_vm, az, el)

if nargin <= 5, az = -40;el = 0;end
% close all;
figure;
%set(h,'Color',[0 0 0]);   
%%xlim([-4 4]);ylim([-4 4]);zlim([-1.2 3]);
view(az,el);axis equal;axis off;

%% generate basis for plotting
nicos = 2;      % number of icosahedron subdivisions used for self-intersection test
[X, C]=surface_mesh.sphere_mesh_gen(nicos);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
L_max = shp_surface.get_L_max(X_o_res_vec(1,:));
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
phys.Y_LK_self  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
%% plot the sequence of shapes
for ix = 1:size(X_o_res_vec,1)
    disp(pos(ix,:));
    X_slave = X_o_res_vec(ix,:)';
    [xclks yclks zclks] = get_xyz_clks(X_slave);
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
    cla;
    plot_shell(X, C,phys.D);
    
    %patch('vertices',X,'faces',phys.C_slave,'FaceColor','red','EdgeColor','k', 'FaceAlpha', 1);
    hold on;
    
%     %% plot infinite plane
%     patch('vertices',[4 4 -1; -4 4 -1;0 -4 -1],'faces',[1 2 3],'FaceColor','b','FaceAlpha', 1);
%     
%     %% plot master
%     s_master = position(s_master,pos(ix,:));
%     phys.X_o_vm = s_master.X_o;
%     [xc yc zc] = get_xyz_clks(phys.X_o_vm);
%     phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for vitteline membrane
%     patch('vertices',phys.X_vm,'faces',phys.C_vm,'FaceColor','green', 'FaceAlpha', 1, 'EdgeColor','none');
%     
    
    hold off
    lighting flat;camlight;
    ylim([0 600]);
    drawnow
    
% % % % % % % % %     [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
% % % % % % % % %     [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
% % % % % % % % %     disp([res_self res_vm]);
%           pause
end
rotate3d on;