function [I] = X2bin(X, dim)
% generate a binary volume I for the object(s) defined by the triangular
% mesh X F. We need the mesh because we might have to interpolate

%% generate the empty image
I = zeros(dim(1),dim(2),dim(3), 'double');

%% set voxels that contain a vertex to one
for ix = 1:size(X,1),
    x = ceil(X(ix,2));
    y = ceil(X(ix,1));
    z = ceil(X(ix,3));
    I(x,y,z) = 1;
end
