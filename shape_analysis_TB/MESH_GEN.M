function [x, y, z, X, C, A, V, v, F_areas, h, H, Eb, da] = mesh_gen(dim, A)
if nargin ==0,dim = 20;A = 12.5664;
elseif nargin ==1, A = 12.5664;end
clc
% [x y z] = sphere(dim-1);
P = partsphere(dim^2);
x = reshape(P(1,:),dim,dim);
y = reshape(P(2,:),dim,dim);
z = reshape(P(3,:),dim,dim);
[t p r] = kk_cart2sph(x,y,z);
[x y z] = kk_sph2cart(t,p,sqrt(A/4/pi));
x = reshape(x,dim^2,1);y = reshape(y,dim^2,1);z = reshape(z,dim^2,1);
X = [x y z];
[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
whos X x
C = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
[A, V, v, F_areas, h, H, Eb, da] = triangulated_props(X, C, 1);
x = X(:,1); y = X(:,2); z = X(:,3);
