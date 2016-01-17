function b = background3D(I)

ds = 50;
b = zeros(size(I));
counter = 0;
for ix = 1:size(I,3)
    b(:,:,ix) = imopen(double(I(:,:,ix)),strel('disk',ds));
end
