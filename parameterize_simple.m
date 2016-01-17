function [X_o] = parameterize_simple(X,F)
%%
%% Author: Khaled Khairy

%%
close all;
[E, L] = edge_info(X,F);
figure;[A, V, v, F_areas_o, h, H, wb, da, N, K, k_g] = triangulated_props(X, F, 1);
F_areas_relative = F_areas_o/sum(F_areas_o);
[t,p,dtline,W] = bijective_map_gen(X, F, L, 1); %%%% Generate Bijective Map
save data_bijective_map_gen_small
%%%% Constrained optimization for uniform parametrization
EqualAreas = 1;
maxiter = 5;
JacCalc = 1;
load data_bijective_map_gen_small;
%[t,p, s, area_res] = uniform_dist_gen(t,p,X, F,L,E, F_areas_relative, maxiter, DMCalc, JacCalc, EqualAreas);save data_uniform_dist_gen_small;
[t,p,A] = newton_steps( t,p,F,300,1,'pombe_01');save data_uniform_dist_gen_small;
%%% Expand in SH
load data_uniform_dist_gen_small;
gdim = 40;L_max = 17;
xL_max = L_max;yL_max = L_max; zL_max = L_max; plotflag = 1;shape_prop = 1;
[A, V, h, v, xclks, yclks, zclks] = sh_projection2(gdim, xL_max, yL_max, zL_max, X(:,1)', X(:,2)', X(:,3)', t', p');view(2);
%     save data_sh_projection_small;
X_o = [xclks(:)' yclks(:)' zclks(:)'];
