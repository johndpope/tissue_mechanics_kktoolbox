function [I, Ix, Iy, Iz, Ii, xvec, yvec, zvec, ivec] = sh4p_plot(GNLK,dimI)
%%% Synthesis of Hyperspherical harmonic parameterization
%%% Generates the 3D intensity volume I of dimensions dim x dim x dim,
%%% given the parametric hyperspherical harmonic coefficients.
%%% Input:
%%%     GNLK: cell array of 4 DNLK structs corresponding to X Y Z and I
%%% Author:
%%%     Khaled Khairy, EMBL-Heidelberg 2007, khairy@embl.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = 1;
% size_limit = 1e8;        % bytes/processor/task
tau = dimI;
dim = tau/2;
pad_vec = ([dim dim dim]);
I = zeros(dimI, dimI, dimI,'single');
I = padarray(I,pad_vec);
av_i = round(size(I,1)/2);av_j = round(size(I,2)/2);av_k = round(size(I,3)/2);
[ix iy iz] = ind2sub(size(I),1:length(I(:)));
x = ix -av_i;y = iy-av_j;z = iz-av_k;  % move coordinate system to center of mass
[t p r] = kk_cart2sph(x,y,z);       %% generate theta, phi and r from the image
indx = sub2ind(size(I),ix, iy, iz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xDNLK = GNLK{1};
yDNLK = GNLK{2};
zDNLK = GNLK{3};
iDNLK = GNLK{4};

[Ix]    = sh4d_plot(xDNLK, tau,[],0)  ;
[Iy]    = sh4d_plot(yDNLK, tau,[],0)  ;
[Iz]    = sh4d_plot(zDNLK, tau,[],0)  ;
[Ii]    = sh4d_plot(iDNLK, tau,[],0)  ;

% the quantities Ix Iy Iz and Ii should completely define the intensity
% volume

%%%%%%%%%%%%%%%%%% Generate the actual intensity volume
fac = 200;
xvec = fac*double(Ix);
yvec = fac*double(Iy);
zvec = fac*double(Iz);
ivec = fac*double(Ii);%rand(length(xvec),1);


% dfig;plot3(xvec,yvec,zvec,'.'); axis equal;     % look at the general geometry of this cloud
% vec = 1:3430;xvec = xvec(vec);yvec = yvec(vec);zvec = zvec(vec);ivec = ivec(vec);

vismethod = 2;


if vismethod==1,  %% use alpha maps to visualize
    %% first we have to interpolate to a regular grid.
    minval = min([min(xvec(:)) min(yvec(:)) min(zvec(:))]);
    maxval = max([max(xvec(:)) max(yvec(:)) max(zvec(:))]);
    d = minval:(maxval-minval)/50:maxval;
    [xi,yi,zi] = ndgrid(d,d,d);
     w = griddatan([xvec(:) yvec(:) zvec(:)], ivec(:),[xi(:) yi(:) zi(:)],'nearest',{''});
    %vi = interp3(xi,yi,zi,ivec,xvec,yvec,zvec);


elseif vismethod == 2,

    fac = 1/100;
    x_pixel = 1/fac;y_pixel = 1/fac;z_pixel =1/fac;
    xdim = dimI;
    ydim = dimI;
    zdim = dimI;
    
    I = zeros(xdim, ydim,zdim);
    ix = round(xvec/x_pixel) + round(xdim/2);iy = round(yvec/y_pixel) + round(ydim/2);iz = round(zvec/z_pixel) + round(zdim/2);
    
    ixval = find(ix<(xdim+1) & ix>0);ix = ix(ixval);iy = iy(ixval);iz = iz(ixval);
    iyval = find(iy<(ydim+1) & iy>0);ix = ix(iyval);iy = iy(iyval);iz = iz(iyval);
    izval = find(iz<(zdim+1) & iz>0);ix = ix(izval);iy = iy(izval);iz = iz(izval);
    izval = find(iz<(zdim+1) & iz>0);ix = ix(izval);iy = iy(izval);iz = iz(izval);
    
    indx = sub2ind(size(I),iy,ix,iz);
    I(indx) = ivec(indx);
    kk_montage(mat2gray(I));impixelinfo;


    %%%%%%%%%%%%% Convolve with PSF if specified/required %%%%%%%
    ffton = 0;
    if ffton
        Sfft = fftn((I),[xdim ydim zdim]);
        C_fft = Sfft.*PSF_fft;%%% convolve with PSF
        I(:,:,:)  = ifftn(C_fft);
        I = fftshift(I);
    end

end













