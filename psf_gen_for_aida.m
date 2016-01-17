function psf = psf_gen_for_aida(dim, s)
%%% psf wrapped around
im =gauss_3d_origin(dim, s);
imfft = fftn(im);
imifft = ifftn(imfft);
psf = ifftshift(imifft);
