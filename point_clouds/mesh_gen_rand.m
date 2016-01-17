function [X, F,x, y, z,  A, V, v, F_areas, h, H, Eb, da] = mesh_gen_rand(dim, A)
if nargin ==0,dim = 30;A = 12.5664;
elseif nargin ==1, A = 12.5664;end
clc

%%%% old approximate method
% % P = partsphere(dim^2);
% % x = reshape(P(1,:),dim,dim);
% % y = reshape(P(2,:),dim,dim);
% % z = reshape(P(3,:),dim,dim);

% % %%%% using subdivisions of icosahedron
if dim>6, dim = 6;end
[X,F]=BuildSphere(dim);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,sqrt(A/4/pi));

X = double([x(:) y(:) z(:)]);


[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
whos X x
C = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
F = C;

[A, V, v, F_areas, h, H, Eb, da] = triangulated_props(X, F, 0);

x = X(:,1); 
y = X(:,2); 
z = X(:,3);
% else
%     error('Too many triangles');
% end
