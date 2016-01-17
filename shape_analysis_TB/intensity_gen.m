function intensity_gen(xclks, yclks, zclks)
%%% given the clks, the basis functions and the PSF, the 3D intensity image is generated

global Iobj Iobj2 Sfft C_fft Y_LK PSF_fft x_pixel y_pixel z_pixel Xmx Xmy Xmz ffton avInt

xdim = size(Iobj, 1);ydim = size(Iobj, 2);zdim = size(Iobj, 3);

gdim = size(Y_LK,1);
Xmx(:) = Y_LK*yclks(:);
Xmy(:) = Y_LK*xclks(:);
Xmz(:) = Y_LK*zclks(:);
%%%%%% Generate the binary volume
Iobj(:) = 0;
% ix = round(Xmx/x_pixel) + round(xdim/2);iy = round(Xmy/y_pixel) + round(ydim/2);iz = round(Xmz/z_pixel) + round(zdim/2);
ix = round(Xmx) + round(xdim/2);
iy = round(Xmy) + round(ydim/2);
iz = round(Xmz) + round(zdim/2);

ixval = find(ix<(xdim+1) & ix>0);ix = ix(ixval);iy = iy(ixval);iz = iz(ixval);
iyval = find(iy<(ydim+1) & iy>0);ix = ix(iyval);iy = iy(iyval);iz = iz(iyval);
izval = find(iz<(zdim+1) & iz>0);ix = ix(izval);iy = iy(izval);iz = iz(izval);
indx = sub2ind(size(Iobj),ix,iy,iz);
Iobj(indx) = 1;



% % Sfft = fftn((Iobj),[xdim ydim zdim]);
% % C_fft = Sfft.*PSF_fft;%%% convolve with PSF
% % Iobj(:,:,:)  = ifftn(C_fft);
% % Iobj = fftshift(Iobj);