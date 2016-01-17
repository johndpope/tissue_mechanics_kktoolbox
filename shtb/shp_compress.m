function clks = shp_compress(I, L_max)
%%% Algorithm for "lossy" image compression using 3D spherical harmonics
%%% for 2D intensity images. The general idea is to map the image
%%% onto the unit sphere, and use the intensity values as the radial
%%% distance. For simple smooth transitions in the image this can be
%%% directly expanded in spherical harmonics (use sh_compress).
%%% For the general case of sharp features in the image, we
%%% should perform the Cartesian parametric expansion. The "shape" is
%%% then interpreted as x y z Cartesian coordinates which are represented
%%% as a series expansion individually. In this case use shp_compress
%%% Usage: clks = sh_compress(I, lmax)
%%%
%%% Limitations: The code is still only for 2D intensity images

verbose = 1;
I = mat2gray(I);
if nargin ==1, L_max = 6;end
clks = [];
%% map image to sphere
n = size(I,1);
gdim = size(I, 1);
[x y z] = sphere(gdim-1);
% surf(x,y,z,I);axis equal;
[t p r] = kk_cart2sph(x,y,z);
r = I;
[x y z] = kk_sph2cart(t,p,r);
%surf(x,y,z);axis equal;

%% do the SH analysis
X = [x(:) y(:) z(:)];
global S U V invS
str = sprintf('SVD_EquiAngularGrid_gdim_%d_lmax_%d.mat', gdim, L_max);
if exist(str)==2, load(str);if verbose, disp('Loading precalculated SVD of SH basis functions for angular grid');end;
    [xclks] = sh_expand(x(:)', L_max, t(:)',p(:)', 0);
else
    [xclks, U, S, V, invS] = sh_expand(x(:)', L_max, t(:)',p(:)', 1);
    save(str,'U', 'S', 'V', 'invS');
end;
[yclks] = sh_expand(y(:)', L_max, t(:)',p(:)', 0);
[zclks] = sh_expand(z(:)', L_max, t(:)',p(:)', 0);
clear S U V invS;

clks = [xclks(:)' yclks(:)' zclks(:)'];

% look at results by decompressing
disp('Decompressing for display');
imshow([I shp_decompress(clks, gdim)]);


% figure;sh_plot_uniform_mesh(clks);
% figure;sh_plot_sphere(clks);



