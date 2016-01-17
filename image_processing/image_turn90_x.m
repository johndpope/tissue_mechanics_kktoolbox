function B2 = image_turn90_x(B, k)
% turn the image k*90 degrees around the spim x axis, i.e. the "row axis" of the 3d
% matrix
if nargin<2, k = 1;end
B = permute(B,[3 2 1]);
B2 = zeros([size(B, 2), size(B, 1), size(B,3)]);
for ix = 1:size(B,3), B2(:,:,ix) = rot90(B(:,:,ix),k);end
B2  = permute(B2,[3 2 1]);