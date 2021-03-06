function I = read_tif_frames(prfx, first, last)
%% usage: I =  load_stack(prfx, first, last); % where prfx is a string 
%% the filename of the stack.
step = 1;
str = sprintf('%s00%d.tif',prfx,first);
I = zeros([size(imread(str)) length(first:step:last)]);
counter = 0;
for ix = first:step:last,
    ix
    counter = counter + 1;
    if ix<10
        str = sprintf('%s00%d.tif',prfx,ix);
    elseif ix<100
        str = sprintf('%s0%d.tif',prfx,ix);
    elseif ix<1000
        str = sprintf('%s%d.tif',prfx,ix);
    end
    I(:,:,counter) = imread(str, 'tiff');
end