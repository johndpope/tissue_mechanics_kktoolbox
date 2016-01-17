function background = background2D(I)

ds = 50;
background = zeros(size(I,1), size(I,2), 1);
background = background + imopen(I,strel('disk',ds));

