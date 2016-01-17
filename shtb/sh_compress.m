function clks = sh_compress(I, L_max)
%%% Algorithm for "lossy" image compression using 3D spherical harmonics is
%%% presented on an example 2D intensity image. The general idea is to map the image
%%% onto the unit sphere, and use the intensity values as the radial
%%% distance. For simple smooth transitions in the image this can be
%%% directly expanded in spherical harmonics.
%%% For the general case of there existing sharp features in the image, we
%%% have to perform the Cartesian parametric expansion. The "shape" is
%%% then interpreted as x y z Cartesian coordinates which are represented
%%% as a series expansion individually.
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
% surf(x,y,z);axis equal;

%% do the SH analysis
X = [x(:) y(:) z(:)];
% [clks] = sh_analysis_LS(X, L_max, 0);
str = sprintf('EquiAngularGrid_gdim_%d_lmax_%d.mat', gdim, L_max);
if exist(str)==2, load(str);if verbose, disp('Loading precalculated SH basis functions for angular grid');end;
    [clks] = sh_expand(r(:)', L_max, t(:)',p(:)', A);
else
    [clks, A] = sh_expand(r(:)', L_max, t(:)',p(:)');
    save(str,'A');
end;


% look at results
disp('Decompressing for display');
imshow([I sh_decompress(clks, gdim)]);

% figure;sh_plot_uniform_mesh(clks);
% figure;sh_plot_sphere(clks);


% P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
% [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
% X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
% cla;patch('Vertices', X, 'Faces', C,'FaceVertexCData',I(:), 'FaceColor', 'interp','EdgeColor','none');
% daspect([1 1 1]);axis off; lighting gouraud;lightangle(64,-42);lightangle(-100,0);view(69,-43);




