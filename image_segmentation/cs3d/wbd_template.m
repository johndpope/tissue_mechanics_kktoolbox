%%%% template for dgvf segmentation
%%%% Requires config_cs3d_t50.mat

disp(' >>>>>>> Diffusion gradient vector field segmentation <<<<<<<<');
load config_cs3d_t50.mat;     % loads parameters good for time 50 and probably for the rest too
%filename = 'I_data_two_ellipsoids.tif';
filename = 'blob_smooth.tif';
I = read_tif_stack(filename);
%write_tif_stack(I,filename);
%% gradient only
[un, vn, wn] = gradient(double(I));
%% call dgvf only
disp(['Generating diffusion gradient vector field ']);
[un, vn, wn] = dgvf_calc(double(I), dgvf_iter, miu, dt, dx, dy, dz, 0); % Note: distributed is set to zero here
disp('Done');
%%
zfac = 1;
dsfac = 1;
L_max = 4;
gdim = 20;
aniso_iter = 0;
field_thresh = 0.0;     % a value of zero takes a lot of time
[segmented, num, V, all_X, all_F, T_final, A, Entropy, Eb] =...
    cs3d(...
    dsfac, zfac,...                                     %downsampling factor and interpolation in z:  keep 1 if coming form cs3d_preprecssing
    distributed,filename,...                            % general configuration
    aniso_iter,delta_t,kappa,option, 1, 1, 1,...        % for anisotropic diffusion
    dgvf_iter, miu, dt, dx, dy, dz,...                  % for dgvf configuration
    gft_niter, field_thresh, ...                        % configure gft
    lbnd,...                                            % configure local adaptive thresholding
    V_max, V_min,...                                    % configure rendering
    laplace_iter,L_max, gdim, newton_iter);