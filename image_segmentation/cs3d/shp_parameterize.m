function [X_o] = shp_parameterize(X, F, L_max, gdim, newton_niter)
dbg = 0;
use_fastrbf_triangulation;
if ~euler_violation,
    
    [t,p,dtline,W] = bijective_map_gen(X, F, L, 0); %%%% Generate Bijective Map
    [t,p] = newton_steps( t,p,F,newton_niter,2,'object',0);%%%% Constrained optimization for uniform parametrization
    [A, V, h, v, xclks, yclks, zclks] = sh_projection2(gdim,L_max, L_max, L_max, X(:,1)', X(:,2)', X(:,3)', t', p');%%% Expand in SHt
    X_o = [xclks(:);yclks(:);zclks(:)];
else
    X_o = [];
    warning('Shape not closed!!');
end
%%%%%%