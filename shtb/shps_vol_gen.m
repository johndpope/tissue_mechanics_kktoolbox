function [I, PSF, voxel] = shps_vol_gen(S, dim)
% Returns a 3D matrix (image volume)
% Usage: [I, PSF, vdim] = shps_vol_gen(T, dim)
% Input:
%       T: usual shape matrix n x m, where n is the number of shapes
%       dim: dimensions of output image, for example [128 128 128]
%           for images where for instance the z-resolution is twice as bad as x and y
%           then the above dim would be [128 128 64]
% Output:
%       I: the 3D intensity image (8-bit)
%       PSF: the point spread function with which the image was formed
%       vdim: voxel dimensions in units of the shape coefficients (for eg. micron)
% Example:
%           [I, PSF, vdim] = shps_vol_gen(T, [64 64 32]);
% Copyright: Khaled Khairy 2009 EMBL
% Also see: write_shps_movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(size(S,2),3)==1, S = S(:,1:end-1);end; % we only want the shape coefficients (not intensity)
I = zeros(dim);
dimfac = 10/100;    % the factor by which the dimension will be extended beyond that of the shps
gdim = 30;
X = [];
%%%%
Lmax = sqrt(size(S,2)/3)-1;
cla;
%%% determine center of mass
disp('Calculating center of mass...');
[Xs,C,Y_LK]=sh_calc(S(1,:), gdim, Lmax);   % get basis functions and connectivity once
for six = 1:size(S,1),
    [Xs]=sh_calc(S(six,:), gdim, Lmax, C, Y_LK);
    X = [X ; Xs];
    Xsmx(:,:,six) = Xs;
    %Csmx(:,:,six) = C;
end
cm = sum(X,1)/length(X);
disp('Done!');
%%% center the shapes
for six = 1:size(S,1),
    [Xs]=sh_calc(S(six,:), gdim, Lmax, C,Y_LK);
    Xs(:,1) = Xs(:,1)-cm(1);          % to center around the center of mass
    Xs(:,2) = Xs(:,2)-cm(2);          % to center around the center of mass
    Xs(:,3) = Xs(:,3)-cm(3);          % to center around the center of mass
    Xsmx(:,:,six) = Xs;               % store all coordinates after translation
    patch('Vertices', Xs, 'Faces', C,'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
end
drawnow;
%%%% Determine voxel size
maxx = max(X(:,1))*dimfac;
maxy = max(X(:,2))*dimfac;
maxz = max(X(:,3))*dimfac;
voxel = [maxx/dim(1) maxy/dim(2) maxz/dim(3)];  % this is the voxel size in space dimensions determined by the CLKs

in = inhull([xi(:) yi(:) zi(:)], X, C);

I = zeros(dim(1),dim(2),dim(3));

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,C, Y_LK]=sh_calc(X_o, gdim, Lmax, C, Y_LK)
X_o = tr(X_o, Lmax);
nc = length(X_o)/3;
xclks = X_o(1:nc);
yclks = X_o(nc+1:2*nc);
zclks = X_o(2*nc+1:end);
if nargin<4,
    P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
    [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
    X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
    Y_LK = get_basis(t',p',gdim,Lmax);
end
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
end

