%%% test bijective mapping
close all;clc;
% load('/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/@surface_mesh/XF_cube.mat');
%load('/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/@surface_mesh/XF_pawn.mat');
%load('/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/@surface_mesh/XF_discocyte.mat');
fn = '/Users/khairyk/VPCZ1_kk_share/mwork/kktoolbox_local/tissue_mechanics_kktoolbox/@surface_mesh/letterE_small.ply';[X,F] = read_ply(fn);
m = surface_mesh(X, F);
m = m.edge_info;
%%% configure
m.newton_niter = 300;
m.newton_step = 0.5;
m.needs_map2sphere = 1;
m.bijective_plot_flag = 1;
m.mapping_plot_flag = 2;
%%%%% Plan:
%%% - calculate local curvature per vertex and also per edge
%%% - generate a data structure with edge weights equal to local curvature
%%% - find the most distant two vertices based on edge weights (ixN and
%%%   ixS)
[G, d, m.ixN, m.ixS] = get_graph(m);
%%
%%% - generate a bijective mapping based on diffusion equation using edge


%%% curvature to penalize/encourage diffusion



mm = map2sphere(m);   % map to sphere

%[t,p,dtline, W] = surface_mesh.bijective_map_gen(X, F, m.L, 1, m.ixN, m.ixS);

% 
% [t, A] = surface_mesh.latitude_calc(m.L, m.ixN, m.ixS); 
% 
% [p, dtline, W, A, b] = surface_mesh.longitude_calc(...
%                        m.X(:,1), m.X(:,2), m.X(:,3),t, A, m.F, m.L,...
%                        m.ixN, m.ixS); 
% 
%                    
% % figure;patch('Vertices',m.X,'Faces',m.F,...
% %     'FaceVertexCData',p,'FaceColor','interp',...
% %     'EdgeColor','k');axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;
% % 
% 
% 
% figure;
% [u, v, w] = kk_sph2cart(t,p,ones(size(p)));
% plot_state(u,v,w,m.F);drawnow

%% parameterize on sphere with spherical harmonics and save
fn = 'current_surface.shp3';
s = shp_surface(mm, 16);
export_ascii(s,fn);
figure;plot_pretty(s, 4);