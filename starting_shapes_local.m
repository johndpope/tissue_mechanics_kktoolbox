function [T, V_vec] = starting_shapes_local(Iw,Idata, L_max, newton_niter, V_max, V_min)
%%% returns a starting set of parameters based on thresholding the data and
%%% estimating CLKs for all found shapes.
%% Important the input matrix I is assumed to be segmented already
% USAGE: [T, A_vec, V_vec] = starting_shapes_local(Isegmented,I, 3, 100,500, 50);
% Iw is the segmented dataset
% Idata is the raw or deconvolved dataset

T = [];
V_vec = [];
disp('Number of segmented objects found:');
disp(max(Iw(:)));
BW = zeros(size(Iw),'uint8');
%for t = min(I(:)):max(I(:))
for t = max(Iw(:)):-1:1%max(Iw(:))-100
    disp('Generating object:'); disp(t);
    %BW = (Iw==t);
    %[L, num] = bwlabeln(BW);
    BW = (Iw==t);BW(BW) = 1; 
    L = bwlabeln(BW);
    STATS = regionprops(L,'Centroid', 'BoundingBox');
    centers = cat(1,STATS.Centroid);
    b = cat(1,STATS.BoundingBox);
    %for ix = 1:num,
        V = sum(sum(sum(L==1)));
        if V<V_max && V>V_min,
            Iw = Iw-(L==1);
            tmp_im = (L==1);
            tmp_im = padarray(tmp_im,[10 10 10]);
            [X, F] = triangulate_gray(squeeze(tmp_im),1, 1,1,.5,size(Iw,4),0,1);%% %% Segmentation and triangulation
            if ~isempty(X)
                newShape = shp_parameterize(X,F,L_max, 30,newton_niter);
                if ~isempty(newShape),
                    [xc, yc, zc] = get_xyz_clks(newShape);
                    px = centers(1,1);py = centers(1,2);pz = centers(1,3); pfac = 3.5449;%1/Y_LK(1);%%
                    xc(1) = px * pfac;yc(1) = py * pfac;zc(1) = pz * pfac;
                    newShape = [xc(:)' yc(:)' zc(:)'];
                    I_roi = Idata(ceil(b(1,2):(b(1,2)+b(1,5))), ceil(b(1,1):(b(1,1)+b(1,4))), ceil(b(1,3):(b(1,3)+b(1,6))));
                    lambda_region = mean(I_roi(:));
                    best_lambda = lambda_region;
                    T = [T;[newShape best_lambda]]; %#ok<AGROW>
                    newShape = [];
                    V_vec = [V_vec V];
                end
            end
        else
            str = sprintf('Object volume = %.2f , \nthis is outside the volume bounds of max %.2f to min %.2f',...
                        V, V_max, V_min);
            disp(str);
        end
    %end
end
disp([V_vec]);
close all;[A, V] = plot_shp_shapes(T(:,1:end-1), 20, 'b');hold off
%close all;[A, V] = plot_shp_T(T, 20, 'b');hold off
save pre_processing_segmentation T V_vec

save bead_draq5_angle0_stvec.mat T;
%analyse_bead_draq5_angle0;

% % % % %%% to scale the whole set back to the original raw data pixel values
% % % % T1 = T;
% % % % sc =  ([1344 1024 1031]./256).*0.32225;
% % % % T2 = T1;
% % % % xa = zeros(1,size(T,1));
% % % % ya = zeros(1,size(T,1));
% % % % za = zeros(1,size(T,1));
% % % % for ix = 1:size(T,1),      %loop over the shapes
% % % %     [xc yc zc] = get_xyz_clks(T(ix,1:end-1));
% % % %     xc = xc.*sc(1);
% % % %     yc = yc.*sc(2);
% % % %     zc = zc.*sc(3);
% % % %     xa(ix) = xc(1);ya(ix) = yc(1);za(ix) = zc(1);
% % % %     T2(ix,1:end-1) = [xc(:)' yc(:)' zc(:)'];
% % % % end
% % % % xm = mean(xa);ym = mean(ya);zm = mean(za);
% % % % 
% % % % for ix = 1:size(T2,1),
% % % %     [xc yc zc] = get_xyz_clks(T(ix,1:end-1));
% % % %     xc(1) = xc(1)-xm;
% % % %     yc(1) = yc(1)-ym;
% % % %     zc(1) = zc(1)-zm;
% % % %     T2(ix,1:end-1) = [xc(:)' yc(:)' zc(:)'];
% % % % end
% % % % 
% % % % % for ix = 1:size(T2,1),
% % % % %     [xc yc zc] = get_xyz_clks(T2(ix,1:end-1));
% % % % %     [t, p, r] = kk_cart2sph(xc(1),yc(1),zc(1));
% % % % %     r = 350;
% % % % %     [u v w] = kk_sph2cart(t,p,r);
% % % % %     xc(1) = u;
% % % % %     yc(1) = v;
% % % % %     zc(1) = w;
% % % % %     T2(ix,1:end-1) = [xc(:)' yc(:)' zc(:)'];
% % % % % end


close all;[A, V] = plot_shp_T(T, 20, 'b');hold on
xlabel('x:micron');ylabel('y:micron');zlabel('z:micron');
axis tight;axis equal;
% % % %%% draw the cytodex sphere also
% % % P = partsphere(60^2);x = P(1,:);y = P(2,:);z = P(3,:);
% % % [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,115);
% % % X = [x(:)-50 y(:)-10 z(:)-40]; [C] = convhulln(X, {'Qt'});
% % % patch('Vertices', X, 'Faces', C,'FaceColor', 'r','FaceAlpha',0.3,'EdgeColor','none');
daspect([1 1 1]);axis on; lighting gouraud;lightangle(64,-42);lightangle(-100,0);view(69,-43);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_o] = shp_parameterize(X,F, L_max,gdim,newton_niter)
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
    filename = 'ddmcmc';
    [t,p,A] = newton_steps( t,p,F,newton_niter,1,filename);save data_uniform_dist_gen;
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

% save data_sh_projection;
% save stvec xclks yclks zclks
% save X_o X_o








