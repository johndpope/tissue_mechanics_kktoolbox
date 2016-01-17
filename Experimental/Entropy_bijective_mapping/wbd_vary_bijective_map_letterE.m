%% mapping letter E
clear all;clc;close all;
fn = 'lettere_quad_option01.off';
[X,F] = read_coff_quad(fn);

%%
Ft = [F(:,[1 2 3]); F(:,[1 4 3])];
m = surface_mesh(X,Ft);
m.laplacian_smooth_iter = 300;
m = laplacian_smooth(m);

[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA, n] = ...
    triangulated_props(m.X, m.F );

% % %% smooth the curvature information
[E,L,face_memb,Ld, Li] = edge_info(m.X,m.F);
Hs = H;
for ix = 1:size(X,1),Hs(ix) = sum(H([L{ix}]))/3;end
for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
% % for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
% % for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
% % for ix = 1:size(X,1),Hs(ix) = sum(Hs([L{ix}]))/3;end
figure;patch('Vertices', X, 'Faces', F,'FaceVertexCData',-Hs, 'FaceColor', 'interp','EdgeColor', 'none','FaceAlpha',1);
%%

m.ixN = get_index(X,[5.8387 0.0866 0.5983]);
m.ixS = get_index(X,[-6.4797 0.3696 -0.1281]);
m = map2sphere(m);
%dfig;[u, v, w] = kk_sph2cart(m.t,m.p,ones(size(m.p)));plot_state(u,v,w,m.F);drawnow