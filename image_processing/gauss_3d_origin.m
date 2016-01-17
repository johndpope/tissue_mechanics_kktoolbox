function G = gauss_3d_origin(dim, s)
% Calculates the 3D Gaussian (for use as a 3D PSF)
% Input:
%       dimx, dimy, dimz: defines the mesh size
%       sx, sy, sz: the standard deviations in the 3 Cartesian directions
%       plotflag: 0 = do not plot slices, 1 = plot slices
%       slice_step: steps taken between outputting slices
%       outputflag:   when ==1 writes to disk tif file with the z-section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
%       G: 3-d matrix of dimensions dim
% The 3D Gaussian is centered around [0 0 0].
% Example:
%
%           G = gauss_3d_origin([128 128 128], [3, 3, 6]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: K. Khairy  ---- July 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimx = dim(1);
dimy = dim(2);
dimz = dim(3);
sx = s(1);
sy = s(2);
sz = s(3);

siz  = ([dimz dimy dimx]-1)/2;
[x, y, z] = meshgrid(-siz(3):siz(3), -siz(2):siz(2), -siz(1):siz(1));
arg = -((x.*x/sx/sx) + (y.*y/sy/sy) + (z.*z/sz/sz));
N = 1/2/pi/sqrt(2 * pi)/sx/sy/sz;
G   = exp(arg);
