function I = get_im_region(filename, w, h, z)
%% load into memory a region of a tif stack file
info = imfinfo(filename);%disp(info);
zdim = length(info);
im = imread(filename,1);

I = zeros([length(h), length(w), length(z)], class(im));
for ix = 1:length(z)
    im = imread(filename,z(ix));
    I(:,:,ix) = im(h,w);
end
