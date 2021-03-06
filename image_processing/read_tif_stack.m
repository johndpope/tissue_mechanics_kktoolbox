function I = read_tif_stack(filename)
%% load a tif stack file
info = imfinfo(filename);%disp(info);
zdim = length(info);
im = imread(filename,1);

if size(im,3)==1,
    if islogical(im),
        I = zeros([size(im,1), size(im,2), zdim], 'uint8');
    else
        I = zeros([size(im,1), size(im,2), zdim], class(im));
    end
    for ix = 1:zdim
        I(:,:,ix) = imread(filename,ix);
    end
    
    
else
    I = zeros([size(im,1), size(im,2), 3], class(im));
    for ix = 1:zdim
        I(:,:,ix) = imread(filename,ix);
    end
end