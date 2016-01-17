function write_tif_frames(I,prfx)
%% usage: save_stacks(I,prfx); % where prfx is a string that will prepend
%% the filename of the stack.

for ix = 1:size(I,3),
    if ix<10
        str = sprintf('%s_00%d.tif',prfx,ix);
    elseif ix<100
        str = sprintf('%s_0%d.tif',prfx,ix);
    elseif ix<1000
        str = sprintf('%s_%d.tif',prfx,ix);
    end
    imwrite(I(:,:,ix),str, 'tiff', 'Compression', 'none');
end