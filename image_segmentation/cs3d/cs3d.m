function [segmented, num, V, all_X, all_F, T_cleaned, A, Entropy, Eb, ILoc] = ...
    cs3d(...
    dsfac, zfac,...                                     % downsample and make isotropic
    distributed,filename,...                          % general configuration
    aniso_iter,delta_t,kappa,option, vxsx, vxsy, vxsz,...  % configure anisotropic diffusion
    dgvf_iter, miu, dt, dx, dy, dz,...                  % configure dgvf
    gft_niter, field_thresh, ...                        % configure gft
    lbnd,...                                            % configure local adaptive thresholding
    V_max, V_min,...                                    % configure surface rendering
    laplace_iter,L_max, gdim, newton_iter)                           % configure SHP fitting              
%%% diffusion gradient vector flow segmentation based on Xu/Prince 1998
%%% the gradient flow tracking based on Li/Wong and thresholding of Otsu
%%% Note the voxel dimensions of I must be isotropic
%%% Usage example:
% % % aniso_iter = 0;             % # of anisotropic diffusion iterations (zero is also ok)
% % % delta_t = 3/44;             % for meaning of remaining anisotropic diffusion parameters please
% % % kappa = 70;                 % type >help anisodiff3D. Also see:
% % % option = 2;                 % Perona and Malik 1990.
% % % voxel_spacing = ones(3,1);
% % % 
% % % dgvf_iter = 200;             % number of diffusion gradient iterations
% % % miu = 5;                    % large values make vector field smoother
% % % dt = 0.005;                 % try to keep dt small (smaller than 0.01 for instance)
% % % dx = 1;                     % don't change this unless you know what your doing
% % % dy = 1;                     % for details see paper: Xu and Prince 1998
% % % dz = 1;
% % % distributed = 1;     % set to zero if Distributed toolbox is not installed or OUT OF MEMORY errors occur
% % % field_thresh = 0.01;       % threshold for points in the vector field to 
% % %                            % consider to propagate (zero is ok but slower)
% % % gft_niter = 40;      % number of voxels estimated between furthest voxel from its sink
% % % lbnd = 5000;        % lowerbound on the size of mosaic region
% % % V_max = 50000;       % Maximum volume to consider for surface rendering
% % % V_min = 200;         % Minimum volume to consider for surface rendering
% % % L_max = 2;          % maximum L order for SHP fitting. Increase for complicated shapes.
% % % gdim  = 30;         % mesh size for L fitting (20->60 is good)
% % % newton_iter = 40;   % number of fitting iterations. Increase for complicated shapes(30->100 is good)
% % % [segmented, mosaic, num, V] = ...
% % %     cs3d(I,distributed,filename,...                     % general configuration
% % %     aniso_iter,delta_t,kappa,option, voxel_spacing,...  % for anisotropic diffusion
% % %     dgvf_iter, miu, dt, dx, dy, dz,...                  % for dgvf configuration
% % %     gft_niter, field_thresh, ...                        % configure gft
% % %     lbnd,...                                            % configure local adaptive thresholding
% % %     V_max, V_min...                                       % configure rendering
% % %     );
%%% Author: Khaled Khairy Copyright EMBL 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% deployment input interpretation
if isdeployed,
    dsfac       = str2num(dsfac);
    zfac        = str2num(zfac);
    distributed = 0;
    filename;
    aniso_iter  = str2num(aniso_iter);
    delta_t     = str2num(delta_t);
    kappa       = str2num(kappa);
    option      = str2num(option);
    vxsx        = str2num(vxsx);
    vxsy        = str2num(vxsy);
    vxsz        = str2num(vxsz);
    dgvf_iter   = str2num(dgvf_iter);
    miu         = str2num(miu);
    dt          = str2num(dt);
    dx          = str2num(dx);
    dy          = str2num(dy);
    dz          = str2num(dz);
    gft_niter   = str2num(gft_niter);
    field_thresh= str2num(field_thresh);
    lbnd        = str2num(lbnd);
    V_max       = str2num(V_max);
    V_min       = str2num(V_min);
    laplace_iter= str2num(laplace_iter);
    L_max       = str2num(L_max);
    gdim        = str2num(gdim);
    newton_iter = str2num(newton_iter);
end
%% %%%%%%%%%%%%%%%%%%%%%%   step -1:    definition of preliminary quantities interpreting input and loading of data
if isdeployed, verbose = 0;else verbose = 0;end
voxel_spacing = [vxsx vxsy vxsz];
if verbose, disp('Working directory is:');disp(pwd);end
V = 0; 
I = read_tif_stack(filename);       %filename = filename(10:end);
if verbose, disp('Size of input image is:');disp(size(I));end
fix(clock)
disp('Starting the segmentation algorithm...');
%% %%%%%%%%%%%%%%%%%%%%%%   step 0:     Calculate isotropic image, interpolate and downsample if necessary   
if~(zfac==1 && dsfac ==1)
if verbose, disp('Resizing image');end
I = image_resize(I,size(I,1)*dsfac, size(I,2)*dsfac, size(I,3)*zfac*dsfac);
end
disp(['Image size is: ' num2str(size(I))]);
disp('Saving configuration.');
%save data_original_I I
save data_configuration dsfac zfac distributed ...
            filename aniso_iter delta_t kappa option ...
            vxsx  vxsy  vxsz voxel_spacing dgvf_iter miu dt dx dy dz...
            gft_niter field_thresh lbnd V_max V_min laplace_iter L_max...
            gdim newton_iter;
%% %%%%%%%%%%%%%%%%%%%%%%   step I:     Anisotropic diffusion
if verbose ==1, figure;kk_montage(mat2gray(I));drawnow;end
if aniso_iter,
    disp(['Performing anisotropic diffusion on ' filename]);
    I = anisodiff3D(I, aniso_iter, delta_t, kappa, option, voxel_spacing);
end
%str_aniso_diffuse = sprintf('data_aniso_diffuse_I_%s.mat', filename);
if verbose, str = sprintf('save data_aniso_diffuse_I_%s.mat I;', filename);eval(str);end
I(:,:,1) = 0;I(:,:,end) = 0;
I = mat2gray(I);
%if verbose, figure;kk_montage(I);drawnow;end
%% %%%%%%%%%%%%%%%%%%%%%%   step II:    dgvf  (diffusion gradient vector field)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Generating diffusion gradient vector field ']);
[un, vn, wn] = dgvf_calc(double(I), 10, 0.5, 0.009, 1, 1, 1, 0); % Note: distributed is set to zero here
%if verbose, str = sprintf('save data_dgvf_%s.mat un vn wn x y z;', filename);eval(str); end%clear Ix Iy Iz;
disp('-------------->  Diffusion gradient vector flow completed!');
%% %%%%%%%%%%%%%%%%%%%%%%   step III:   gft (gradient flow tracking)
[L num xs ys zs]  = gft_calc(un, vn, wn,field_thresh, gft_niter, distributed); %#ok<NASGU,NASGU>
if verbose, str = sprintf('save data_gft_%s.mat L num;', filename);eval(str);end
%if verbose, figure;kk_montage(mat2gray(mat2gray(mosaic) + mat2gray(I)));end
clear un vn wn x y z
disp('--------------> Gradient flow tracking completed!');
%% %%%%%%%%%%%%%%%%%%%%%   step III and 1/2: associate voxels with sinks
mosaic = assign_sinks_to_voxels(L,xs,ys,zs);
%% %%%%%%%%%%%%%%%%%%%%%%   step IV:    local adaptive thresholding
[segmented mosaic] = local_adaptive_thresholding_03(I, mosaic, num, lbnd);
num_ren = max(segmented(:));num = num_ren;

str = sprintf('save data_lat_%s.mat segmented;', filename);eval(str);

if verbose==2 close all hidden;end
if verbose, disp('--------------> Local adaptive thresholding completed!'); end
%% %%%%%%%%%%%%%%%%%%%%%%   step V:     Render surfaces
[all_X, all_F, num, V] = surface_render(segmented,V_max, V_min); 
str = sprintf('save data_render_%s.mat all_X all_F num V;', filename);eval(str);

if verbose, disp('-------------- Surface triangulation rendering completed !!');end

%%% plotting and printing rendered
if verbose,
    str_render = sprintf('data_render_%s.mat',filename);
    load(str_render);
    % h = figure('visible','off');
    h = figure;plot_all_X(all_X, all_F);
    fileout = sprintf('intermediate_result_rendering_object_surfaces_%s.fig',filename);
    saveas(h,fileout);
    print -dtiff -r600 intermediate_result_rendering_object_surfaces.tif;
    if verbose==2 close(h);end
end
%% %%%%%%%%%%%%%%%%%%%%%%   step VI:    Process surfaces (repair and clean meshes) --- optional
% % if~isdeployed
% %     global ISO2MESH_SESSION; ISO2MESH_SESSION= filename;
% %     global ISO2MESH_TEMP;ISO2MESH_TEMP = './';
% %     tic;[all_X, all_F] = clean_all_XF(all_X, all_F, laplace_iter);
% %     str_render_clean = sprintf('data_render_clean_%s.mat', filename);
% %     str = sprintf('save data_render_clean_%s.mat all_X all_F;', filename);eval(str);
% %     toc;disp('-------------- Surface triangulation cleaning completed !!');
% %     num = length(all_X);
% %     %%% plotting and printing rendered and cleaned
% %     if verbose
% %         str_render_clean = sprintf('data_render_clean_%s.mat', filename);
% %         load(str_render_clean);
% %         % h = figure('visible','off');
% %         h = figure;plot_all_X(all_X, all_F);axis equal;lighting phong; light;
% %         fileout = sprintf('intermediate_result_rendering_cleaned_object_surfaces_%s.fig',filename);
% %         saveas(h,fileout);
% %         print -dtiff -r600 intermediate_result_rendering_object_surfaces.tif;
% %         if verbose==2 close(h);end
% %     end
% % end
%% %%%%%%%%%%%%%%%%%%%%%%   step VII:   Spherical Harmonics Parameterization
% % load data_configuration
% % str = sprintf('load data_render_%s.mat;', filename);eval(str);

[T, T_c, V_vec, A_vec, wb_vec, k_g_vec, h_vec]  = ...
    shp_gen(all_X, all_F, L_max, gdim, newton_iter, V_max, V_min);

[T_cleaned A V Eb] = clean_T(T_c, A_vec, V_vec, wb_vec);
 
str = sprintf('save data_shp_%s.mat T T_cleaned A V Eb V_vec A_vec wb_vec k_g_vec h_vec;',filename);eval(str);


if verbose, disp('-------------- SHP fitting completed !!');end
if verbose %%% plotting and printing SHP surfaces
    str_shp = sprintf('data_shp_%s.mat',filename);
    load(str_shp);
    h = figure('visible','on');
    [A, V, Y_LK, C] = plot_shp_shapes(T_cleaned(:,1:end-1),20,'random');
    title('Random color coding.');
    saveas(h,'intermediate_result_SHP_energy.fig');
    print -dtiff -r600 intermediate_results_SHP_energy.tif;
    if verbose == 2;close(h);end
end

% % 
%% %%%%%%%%%%%%%%%%%%%%%%   step VIII:  Object region entropy scoring
% str = sprintf('load data_shp_%s.mat;',filename);eval(str);
% str = sprintf('load data_aniso_diffuse_I_%s.mat;', filename);eval(str);

[T_final Entropy delidx ILoc]  = score_T(I, T_cleaned);
A(find(delidx==1)) = [];V(find(delidx==1)) = [];Eb(find(delidx==1)) = [];

num = size(T_final,1);
str = sprintf('save data_score_%s.mat T_final Entropy A V Eb delidx;',filename);eval(str);

if verbose,
disp('-------------- Object scoring completed !!');
disp('----Main segmentation and object parameterization completed ----');
end
