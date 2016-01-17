function Iout = image_3D_fuse(Ic1, Ic2, wavelet_type, level)
%% fuse 3D dataset by looping over the 2D frames and doing the fusion in wavelet space
if nargin == 2,
wavelet_type = 'db4';
level = 5;
end

tic;
Iout = zeros(size(Ic1), class(Ic1));
for ix = 1:size(Ic1,3)
    %disp(['Frame: ' num2str(ix) ' of ' num2str(size(Ic1,3))]);
    im1 = Ic1(:,:,ix);
    im2 = Ic2(:,:,ix);
	Iout(:,:,ix) = uint16(wfusimg(im1,im2,wavelet_type,level,'mean','max'));
end