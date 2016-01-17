function [sub_filenames, zfac, dsfac]= cs3d_preprocessing(filename,MEM_MAX, PADDING)
% partitions the image (filename) into 3D regions that are size-wise small
% enough that a node with 2GB of RAM can run cs3d on them individually
% In the future: preprocessing will also include deconvolution of the image
% partitions.
% Input:
%       filename: name of the 3D tif stack (background subtract and medfilt
%                  in fiji before entering here)
%       camera_chip_size: in micron (e.g. 6.45)
%       objective_magnification: e.g. 10
%       z-spacing: in micron
dsfac = 1;zfac = 1;
%% load the image info
info = imfinfo(filename);zdim = length(info);
imx = info.Width;
imy = info.Height;
imz = length(info);
disp('Dimensions of input image:');disp([imx imy imz]);
%% Esimate the size of the resized isotropic stack (before resizing it)
mem_size = dsfac*(prod([imx imy imz].*[1 1 zfac])*8/1e6);   % the amount of RAM needed to keep

%% partion the stack into individual images by cropping xy windows along z
% (i.e. z always stays the same)
level = 1;
sub_im_size = mem_size;
nx = ceil(imx/(level));          % x-dimension if subdivided at this level
ny = ceil(imy/(level));          % y-dimension at this level
pad = (nx+ny)*2*PADDING/1e6;       % assuming a general region in the middle that gets padded in Mega Bytes
sub_im_size = dsfac*(prod([nx ny imz].*[1 1 zfac])*8/1e6 + pad); % in Mega Bytes

while (sub_im_size>MEM_MAX),        % if it is too large, determine subdivision level (with overlap and recalculate size)
    nx = ceil(imx/(level));          % x-dimension if subdivided at this level
    ny = ceil(imy/(level));          % y-dimension at this level
    pad = (nx+ny)*2*PADDING/1e6;       % assuming a general region in the middle that gets padded in Mega Bytes
    sub_im_size = dsfac*(prod([nx ny imz].*[1 1 zfac])*8/1e6 + pad); % in Mega Bytes
    disp(['Considering level: ' num2str(level) ' with sub-image size (MB): '  num2str(sub_im_size)]);
    level = level+1;
end
disp(['We will need to subdivide to level: ' num2str(level-1)]);
level = level-1;
%% now that we determined the level, let us subdivide, save the images (giving them sequential names) and
% flush memory of the original (large) image. We have to keep track of
% information to reassemble the segmented objects spatially
depth_range = 1:imz;
sub_filenames = {};
counter = 0;
if level==0, level = 1;end  % cluge alarm
for ix = 0:level-1     % iterate over the subdivisions
    height_range = (ny*(ix) + 1 - PADDING):(ny*(ix+1))+ PADDING;         % the y (height) range of voxels to take into this sub-image
    height_range(height_range<1) = [];height_range(height_range>imy) = [];
    for jx = 0:level-1
        counter = counter+1;
        sub_fn = sprintf('%s_%s.tif', filename, get_suffix(3,counter));      % the name of the sub-image file
        disp(['Processing sub image: ' sub_fn]);
        width_range = (nx*(jx) + 1 - PADDING):(nx*(jx+1)+ PADDING);         % the x range of voxels to take into this sub-image
        width_range(width_range<1) = [];width_range(width_range>imx) = [];
        im = mat2gray(get_im_region(filename,width_range, height_range, depth_range));
        im = image_resize(im,size(im,1)*dsfac, size(im,2)*dsfac, size(im,3)*zfac*dsfac);   %resize the images to make them isotropic
        %im = padarray(im,[5 5 5]);
        write_tif_stack(im,sub_fn);
        sub_filenames{counter} = sub_fn;
        data_sub_fn = sprintf('save data_%s_%s.mat width_range height_range imx imy imz zfac dsfac', filename, get_suffix(3,counter));eval(data_sub_fn);
    end
end

%% just in case sub_filenames screws up generate the list here
% % filename = 'fp_probe1_angle0_ch0.tif';%%%%%% Imaging conditions and microscope configuration
% % counter = 0;
% % level = 10;%
% % zfac = 1;
% % for ix = 0:level-1     % iterate over the subdivisions
% %     for jx = 0:level-1
% %         counter = counter+1;
% %         sub_fn = sprintf('%s_%s.tif', filename, get_suffix(3,counter));      % the name of the sub-image file
% %         sub_filenames{counter} = sub_fn;
% %     end
% % end
% % save sub_filenames sub_filenames zfac dsfac