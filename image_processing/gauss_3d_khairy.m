function G = gauss_3d_khairy(dimx, dimy, dimz,sx, sy, sz,plotflag, slice_step, outputflag)
% Calculates the 3D Gaussian (for use as a 3D PSF)
% Input:
%       dimx, dimy, dimz: defines the mesh size
%       sx, sy, sz: the standard deviations in the 3 Cartesian directions
%       plotflag: 0 = do not plot slices, 1 = plot slices
%       slice_step: steps taken between outputting slices
%       outputflag:   when ==1 writes to disk tif file with the z-section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
%       G: 3-d matrix of dimensions dimx, dimy, dimz
% The 3D Gaussian is centered around [0 0 0].
% Example:
%
%           G = gauss_3d_khairy(128, 128, 128, 30, 30, 90, 1, 5, 1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: K. Khairy  ---- July 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

siz  = ([dimz dimy dimx]-1)/2;
[x, y, z] = meshgrid(-siz(3):siz(3), -siz(2):siz(2), -siz(1):siz(1));
arg = -((x.*x/sx/sx) + (y.*y/sy/sy) + (z.*z/sz/sz));
N = 1/2/pi/sqrt(2 * pi)/sx/sy/sz;
G   = exp(arg);
% G(G<eps*max(G(:))) = 0;
% G   = N.* exp(arg);
% G = G./max(max(max(G)));
% Plotting starts here
if plotflag
    h = slice(x,y,z,G,[ 0 ],[],[0]);
    alpha('color');
    set(h,'EdgeColor','none','FaceColor','interp',...
        'FaceAlpha','interp');
    alphamap('rampdown');
    alphamap('increase',.1);
    colormap(hsv);
    daspect([1 1 1]);
end
%%%%%%%%%%%%%%% Produce slices if necessary
if outputflag,
    counter = 0;
    for ix = 1:slice_step:size(G,3)
        counter = counter + 1;
        if counter < 10
            filename = sprintf('gaussian_psf_slice_00%d.tif',counter);
        elseif counter <100
            filename = sprintf('gaussian_psf_slice_0%d.tif',counter);
        elseif counter <1000
            filename = sprintf('gaussian_psf_slice_%d.tif',counter);
        end
        imwrite(G(:,:,ix),filename);
    end
end