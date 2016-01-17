function [s m] = read_image2shp(L_max)
% reads the _para file (in vtk format) as it comes as output from image2shp
% and converts it into theta phi and plots them
% example:
% cd C:\KK_share\c_code\vs_projects\image2shp\x64\Debug
% fn = '_para'; 
% [t p] = read_image2shp(fn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0, L_max = 16;end
%% read the positions of the parametric points on the sphere
fn = '_para';
[F X] = readVTK(fn);
%plot_mesh(X,F);
[t p r] = kk_cart2sph(X(:,1), X(:,2), X(:,3));
[x y z] = kk_sph2cart(t, p, 1);
X = [x(:) y(:) z(:)];
figure;
plot_mesh(X,F);title('position of parametric points');

%% now read the actual surface fields x y z
fn = '_surf';
[F X] = readVTK(fn);
figure;
plot_mesh(X,F);title('surface from voxels');

%% now generate the spherical harmonics parameterization object
m = surface_mesh(X,F);
m.t = t;
m.p = p;
m.needs_map2sphere = 0;
s = shp_surface(L_max,m);