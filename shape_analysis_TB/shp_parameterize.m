function [X_o] = shp_parameterize(X,F, L_max,gdim,newton_niter)
L = {};
use_fastrbf_triangulation;save data_triangulate_gray;
if ~euler_violation,
    %%%%%%%%%
    %%%% Generate Bijective Map
    load data_triangulate_gray;[t,p,dtline,W] = bijective_map_gen(X, F, L, 0); save data_bijective_map_gen
    %%%% Calculate shape properties based on triangulation (optional)
    [A, V, v, F_areas_o] = triangulated_props(X, F, 0);save data_triangulated_props

    %%%% Constrained optimization for uniform parametrization
    load data_triangulated_props;load data_bijective_map_gen;
    % [t,p] = uniform_dist_gen(t,p,X, F, L, E,F_areas_o/sum(F_areas_o),50, 0, 1, 1);save data_uniform_dist_gen;
    filename = 'untitled';
    [t,p] = newton_steps( t,p,F,newton_niter,1,filename);save data_uniform_dist_gen;
    %%% Expand in SHt
    load data_uniform_dist_gen;
    % gdim = 50;L_max = 6;

    plotflag = 0;shape_prop = 1;
    [A, V, h, v, xclks, yclks, zclks, X, C,t,p, Y_LK] = ...
        sh_projection2(gdim,L_max, L_max, L_max, X(:,1)', X(:,2)', X(:,3)', t', p', 0);
    X_o = [xclks(:);yclks(:);zclks(:)];
else
    X_o = [];
end