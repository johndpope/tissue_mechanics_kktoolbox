function I = shp_decompress(clks, gdim)
%%% Decompress images encoded by shp_compress
%%% Usage: I = shp_compress(xc,yc,zc, gdim)
%%% 
%%% Limitations: The code is still only for 2D intensity images
I = [];
[xc yc zc] = get_xyz_clks(clks);
L_max = round((sqrt(length(xc))-1));
if nargin ==1, gdim = 128;end
%% generate the original mesh on the sphere
[x y z] = sphere(gdim-1);
[t p r] = kk_cart2sph(x,y,z);
Y_LK = get_basis(t(:),p(:),gdim,L_max);

Ix = Y_LK(:,1:length(xc))*xc(:);
Iy = Y_LK(:,1:length(yc))*yc(:);
Iz = Y_LK(:,1:length(zc))*zc(:);

[t p I] = kk_cart2sph(Ix,Iy,Iz);
I = mat2gray(reshape(I, gdim, gdim));
if nargin == 1, imshow(I);end






