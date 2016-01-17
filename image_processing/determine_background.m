function bout = determine_background(I)

nbins = 100;
ds = 50;
b = zeros(size(I,1), size(I,2));
counter = 0;
for ix = 1:10:size(I,3)
    counter = counter + 1;
    b = b + imopen(double(I(:,:,ix)),strel('disk',ds));
end
bout = mean(b(:)/counter);
% % disp(b);
% % 
% % [n,xout] = hist(I(:), nbins);
% % disp(xout(find(n==max(n))));