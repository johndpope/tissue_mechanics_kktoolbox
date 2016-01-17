function I = shift_correction_2D(I)
%% corrects for shifts in a 3D stack of images Is
for ix = 2: size(I,3),  % loop over z
    im1 = I(:,:,ix-1);
    im2 = (I(:,:,ix));
    [output Greg] = dftregistration(fft2(im1),fft2(im2),100);
    im3 = (abs(ifft2(Greg)));
    I(:,:,ix) = im3;
end