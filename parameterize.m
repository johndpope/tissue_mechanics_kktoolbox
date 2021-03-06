function [X_o] = parameterize(X,F)
%% use curvature-based mapping onto the parametric sphere for the
%% triangular mesh X F
%% Author: Khaled Khairy
%% Note: This must still be tested
%%

%%
close all;
[E, L] = edge_info(X,F);
figure;[A, V, v, F_areas_o, h, H, wb, da, N, K, k_g] = triangulated_props(X, F, 1);
% %% curvature-based mesh simplification
% %%%we need to identify the vertices of low local mean curvature
% thresh = 1e-9;
% % for ix = 1:10,
% ix_remove = find(abs(K)<= thresh);% generate the list of indices into X to be removed
% disp(length(ix_remove));
% % and remove them
% [X, F] = mesh_simplify(X,F,ix_remove);
% % then recalculate the triangular mesh properties to continue
% %use_fastrbf_triangulation;
% [E, L] = edge_info(X,F);
% figure;[A, V, v, F_areas_o, h, H, wb, da, N, K, k_g] = triangulated_props(X, F, 3);
% %%
%% assign an area factor according to the local Gaussian curvature
cbm_area_fac = sum(abs(K(F)),2)/3/length(F);
cbm_area_fac = cbm_area_fac/norm(cbm_area_fac);
F_areas_relative = cbm_area_fac/sum(cbm_area_fac)*4*pi;
cla; patch('Vertices', X, 'Faces', F,'FaceVertexCData',F_areas_relative, 'FaceColor', 'flat','FaceAlpha',1); axis tight; axis equal;
%%%%%%
% F_areas_relative = F_areas_o/sum(F_areas_o);
 [t,p,dtline,W] = bijective_map_gen(X, F, L, 1,573,272); %%%% Generate Bijective Map
 save data_bijective_map_gen_small
%%%% Constrained optimization for uniform parametrization
EqualAreas = 0;
maxiter = 5;
JacCalc = 1;
for dm =  [-3.8]
    load data_bijective_map_gen_small;
%     DMCalc = 9e-3;JacCalc = 0;
    DMCalc = 1*10^(dm);
%     DMCalc = 0;
    %[t,p, s, area_res] = uniform_dist_gen(t,p,X, F,L,E, F_areas_relative, maxiter, DMCalc, JacCalc, EqualAreas);save data_uniform_dist_gen_small;
     
    [t,p,A] = newton_steps( t,p,F,300,0,'letterE_small');save data_uniform_dist_gen_small;
      [x y z] = kk_sph2cart(t,p, 1);plot_mesh([x(:) y(:) z(:)],F);   
    %%% Expand in SH
    
    load data_uniform_dist_gen_small;
    gdim = 40;L_max = 10;
    xL_max = L_max;yL_max = L_max; zL_max = L_max; plotflag = 1;shape_prop = 1;
     [A, V, h, v, xclks, yclks, zclks] = sh_projection2(gdim, xL_max, yL_max, zL_max, X(:,1)', X(:,2)', X(:,3)', t', p');view(2);
%     save data_sh_projection_small;
    X_o = [xclks(:)' yclks(:)' zclks(:)'];
end
