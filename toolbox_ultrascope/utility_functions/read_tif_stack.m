function I = read_tif_stack(filename)
%% load a tif stack file
info = imfinfo(filename);%disp(info);
zdim = length(info);
im = imread(filename,1);

if size(im,3)==1,
    I = zeros([size(im,1), size(im,2), zdim], class(im));
    for ix = 1:zdim
        I(:,:,ix) = imread(filename,ix);
    end
    
    
else
    I = zeros([size(im,1), size(im,2), 3], class(im));
    for ix = 1:zdim
        I(:,:,ix) = imread(filename,ix);
    end
end