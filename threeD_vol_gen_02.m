function [I, PSF] = threeD_vol_gen_02(X_o, opts)%x_pixel, y_pixel, z_pixel, flag)
%%%%% generate a 3D volume for the shape described by X_o
%%%%% returns the 3D matix of intensities at the boundaries
global Iobj Iobj2 Sfft C_fft Y_LK PSF_fft x_pixel y_pixel z_pixel Xmx Xmy Xmz
global RAW ffton plotflag parm_ix full_parms avInt I_display bkg dim2 dim3

[Area,V,v_red,t,p,X_coord]=plot_sh(X_o,20);
[xclks, yclks, zclks] = get_xyz_clks(X_o);
if nargin == 1,
    opts.dim = [64 64 64];
    vox = max([max(abs(xclks)) max(abs(yclks)) max(abs(zclks))])/32; 
    opts.voxel = [vox vox vox];
%     opts.voxel = [max(abs(xclks)) max(abs(yclks)) max(abs(zclks))]/32;
    opts.conv = 0;
    opts.ill_decrease = 0;
    opts.noise = 0;
    opts.PSF_parms = [1 1 1];
end

dim     = opts.dim;
x_pixel = opts.voxel(1);
y_pixel = opts.voxel(2);
z_pixel = opts.voxel(3);
flag    = opts.conv;
psf_parms = opts.PSF_parms;


Iobj = zeros(dim(1),dim(2),dim(3));
[xclks yclks zclks nc] = get_xyz_clks(X_o);

%% make theoretical PSF
PSF = gauss_3d_khairy(dim(1),dim(2),dim(3),...
    psf_parms(1),psf_parms(2),psf_parms(3),0,5,0);%
PSF_fft = fftn(PSF,[dim(1) dim(2) dim(3)]);

%%%%%%%%%%%%%%%% configure %%%%%%%%%%%%%%%%%%%%
L_max = sqrt(length(xclks))-1;

full_parms = [xclks(:);yclks(:);zclks(:)];
nparms = (L_max+1)^2;
nc = length(xclks);
parm_ix = [1:nparms nc+1:nc+nparms 2*nc+1:2*nc+nparms];
%parm_ix  = 1:length(full_parms);%[3 13];
X_o = full_parms(parm_ix);


%% make surface
finedim = 200;
[Area,V,v_red,t,p,Xcoord,F,Y_LK] = plot_sh(xclks, yclks, zclks, finedim);axis equal
%%%%% prepare initial quantities
Xmx = zeros(length(Xcoord),1);Xmy = Xmx; Xmz = Xmx;
Iobj = zeros(size(PSF_fft));Iobj2 = Iobj;C_fft = zeros(size(PSF_fft));Sfft = zeros(size(PSF_fft));
I_display = zeros(size(PSF_fft,1), size(PSF_fft,2),3,size(PSF_fft,3));

%% make the volume simulatioin
clks = full_parms;
clks(parm_ix) = X_o;
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
ffton = flag;
intensity_gen(xclks, yclks, zclks); % the actual volume generation
I = Iobj;
I = shiftdim(I,2);      % make y to z : so we have (z,x,y)
if opts.ill_decrease    %%% calculate the decrease in illumination along y
%     damp = damp(ones(128,1),:);
%     damp = damp(:,:,ones(128,1));
%     damp = shiftdim(damp,1);

    for slice = 1:size(I,3),
        y_slice = I(:,:,slice); % get the Y-slice
        for p = 64:size(y_slice,2),
            y_slice(p,:) = process_profile(y_slice(p,:));%y_slice(p,:);%
%             y_slice(p,:) = y_slice(p,:);%
        end
        I(:,:,slice) = imrotate(y_slice, 180);%figure(1);imshow(mat2gray(y_slice));drawnow;
    end
%      I = shiftdim(I,1);      % 
    
end


minI = min(I(:));
maxI = max(I(:));
if opts.noise
    %% make RAW from the volume simulation by adding noise and calculating PSF
    Iobj = Iobj + maxI*opts.poisson_bkgrd;
    for ix = 1:size(Iobj,3),
        Iobj(:,:,ix) = Iobj(:,:,ix).*poissrnd(opts.poisson_lambda,size(Iobj,1),size(Iobj,2));
    end;
end
%% use imfill if needed
% Iobj = imfill(Iobj);
%%%


save Iobj_final Iobj;

I = mat2gray(I,[minI (maxI+maxI*0.2)]);
I = gray2ind(I,256);
kk_montage(mat2gray(I));

