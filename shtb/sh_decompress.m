function I = sh_decompress(c, gdim)
%%% Decompress images encoded by sh_compress
%%% Usage: I = sh_compress(clks, gdim)
%%% 
%%% Limitations: The code is still only for 2D intensity images
I = [];
L_max = (sqrt(length(c))-1);
if nargin ==1, gdim = 128;end
%% generate the original mesh on the sphere
[x y z] = sphere(gdim-1);
[t p r] = kk_cart2sph(x,y,z);
Y_LK = get_basis(t(:),p(:),gdim,L_max);
I = Y_LK(:,1:length(c))*c(:);
I = mat2gray(reshape(I, gdim, gdim));






