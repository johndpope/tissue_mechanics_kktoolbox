function crop_frames(directory, prfx, extension, nframes, x, y)
here = pwd;
cd(directory);
for ix = 1:nframes,
    if ix<10
        str = sprintf('%s00%d.%s',prfx, ix,extension);
    elseif ix<100
        str = sprintf('%s0%d.%s',prfx,ix, extension);
    elseif ix<1000
        str = sprintf('%s%d.%s',prfx, ix, extension);
    end
    im = imread(str, extension);
    imcropped = im(y(1):y(2), x(1):x(2), :);
    imwrite(imcropped, str, extension);
end
cd pwd
end