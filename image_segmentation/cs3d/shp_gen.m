function [T, T_cleaned, V_vec, A_vec, wb_vec, k_g_vec, h_vec] = shp_gen(...
    all_X, all_F, L_max, gdim, newton_niter, V_max, V_min)
%% Author: Khaled Khairy. Copyright EMBL-Heidelberg 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1;
T = [];
T_cleaned = [];
V_vec = [];
A_vec = [];
wb_vec = [];
k_g_vec = [];
h_vec = [];
if verbose,disp(['Number of segmented objects found: ' num2str(length(all_X))]);end
T = zeros(length(size(all_X,2):-1:1), 3*(L_max+1).^2);
parfor t = 1:size(all_X,2)
    X = all_X{t};
    F = all_F{t};
    if verbose,disp(['Generating object : ' num2str(t)]);end
    centers = sum(X,1)./size(X, 1);
    newShape = shp_parameterize(X,F,L_max, gdim,newton_niter);
    if ~isempty(newShape),
        [xc, yc, zc] = get_xyz_clks(newShape);
        px = centers(1,1);py = centers(1,2);pz = centers(1,3); pfac = 3.5449;%1/Y_LK(1);%%
        xc(1) = px * pfac;yc(1) = py * pfac;zc(1) = pz * pfac;
        newShape = [xc(:)' yc(:)' zc(:)'];
        T(t,:) = [newShape]; %#ok<AGROW>
    end
end
%[A, V_vec] = plot_shp_shapes(T(:,1:end-1), 20, 'b');hold off

parfor ix = 1:size(T,1),      % loop over the shapes
    if verbose, disp(['Calculating shape properties for object : ' num2str(ix)]);end
    [X,F]=shp_get_coord(T(ix,:), 20);
    [A, V, v, F_areas, h, H, wb, da, N, K, k_g] = triangulated_props(X, F, 0);
    A_vec(ix) = A;
    V_vec(ix) = V;
    v_vec(ix) = v;
    h_vec(ix) = h;
    wb_vec(ix) = wb;
    k_g_vec(ix) = k_g;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



