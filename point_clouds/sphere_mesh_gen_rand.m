function [X, C] = sphere_mesh_gen_rand(dim, A)
if nargin ==0,dim = 120;A = 12.5664;
elseif nargin ==1, A = 12.5664;end

%% old approximate method
P = partsphere(dim^2);
x = reshape(P(1,:),dim,dim);
y = reshape(P(2,:),dim,dim);
z = reshape(P(3,:),dim,dim);
X = double([x(:) y(:) z(:)]);

[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
C = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);


