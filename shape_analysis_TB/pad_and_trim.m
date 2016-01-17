function [RAW2, PSF2,xdim, ydim, zdim, x_pixel, y_pixel, z_pixel] = pad_and_trim(RAW, PSF, x_pixel, y_pixel, z_pixel, fac)
% pads and trims and changes the resolution of RAW and PSF according to the factor fac
% Always returns powers of 2 dimensions suitable for the fft later on


for ix = 1:size(PSF,4), 
    PSF2(:,:,1,ix) = imresize(PSF(:,:,1,ix),fac);
end
for ix = 1:size(RAW,4),     
    RAW2(:,:,1,ix) = imresize(RAW(:,:,1,ix),fac);
end
RAW = RAW2;PSF = PSF2;
x_pixel = x_pixel * 1/fac;y_pixel = y_pixel* 1/fac;
%%% Now that we resized, we need to pad up to the next higher multiple of 2

xdim = max([size(RAW,1) size(PSF,1)]);  %get the xdim that is larger (we assume equal xdim and ydim)
zdim = max([size(RAW,4) size(PSF,4)]);  % get the zdim that is larger
mult2 = [16 32 64 128 256 512]; % reasonable multiples of 2 that can be considered
oneup = find(xdim<mult2);xdim = mult2(oneup(1));    % take the next one up from mult2
oneup = find(zdim<mult2);zdim = mult2(oneup(1));    % take the next one up from mult2
ydim = xdim;

%let's pad the PSF
xdiff = xdim-size(PSF,1);ydiff = ydim-size(PSF,2);zdiff = zdim-size(PSF,4);
if mod(xdiff,2)==0, PSF = padarray(PSF,[xdiff/2 0 0 0]);else PSF = padarray(PSF, [(xdiff-1)/2 0 0 0]);PSF = padarray(PSF, [1 0 0 0], 0, 'post');end
if mod(ydiff,2)==0, PSF = padarray(PSF,[0 ydiff/2 0 0]);else PSF = padarray(PSF, [0 (ydiff-1)/2 0 0]);PSF = padarray(PSF, [0 1 0 0], 0, 'post');end
if mod(zdiff,2)==0, PSF = padarray(PSF,[0 0 0 zdiff/2]);else PSF = padarray(PSF, [0 0 0 (zdiff-1)/2]);PSF = padarray(PSF, [0 0 0 1], 0, 'post');end

%let's pad RAW
xdiff = xdim-size(RAW,1);ydiff = ydim-size(RAW,2);zdiff = zdim-size(RAW,4);
if mod(xdiff,2)==0, RAW = padarray(RAW,[xdiff/2 0 0 0]);else RAW = padarray(RAW, [(xdiff-1)/2 0 0 0]);RAW = padarray(RAW, [1 0 0 0], 0, 'post');end
if mod(ydiff,2)==0, RAW = padarray(RAW,[0 ydiff/2 0 0]);else RAW = padarray(RAW, [0 (ydiff-1)/2 0 0]);RAW = padarray(RAW, [0 1 0 0], 0, 'post');end
if mod(zdiff,2)==0, RAW = padarray(RAW,[0 0 0 zdiff/2]);else RAW = padarray(RAW, [0 0 0 (zdiff-1)/2]);RAW = padarray(RAW, [0 0 0 1], 0, 'post');end

RAW2 = RAW;PSF2 = PSF;
disp(size(PSF));disp(size(RAW));disp([x_pixel y_pixel z_pixel]);