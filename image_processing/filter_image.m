function filteredStack = filter_image(stack,filterMode, r)
%%% filter the stack with filter defined in filterMode, and with range defined in r
%%% legal filter modes are 0 (median), 1 (mean), and 2 (Gaussian);
%%% example: I = filter_image(I,2, 50);
%%% Uses mex file for Gaussian filtering
xSize = size(stack, 1);
ySize = size(stack, 2);
zSize = size(stack, 3);
stack = uint16(mat2gray(stack) *2^16);
switch filterMode
    case 0
        filteredStack = stack - medfilt3(stack, r);
    case 1
        filteredStack = stack - uint16(smooth3(stack, 'box', r));
    case 2
        gaussStack = zeros(xSize, ySize, zSize, 'uint16');
        
        splitting = 5;
        splittingMargin = r + 1;
        
        for i = 1:splitting
            xSlabStart = max(1, round((i - 1) * xSize / splitting + 1 - splittingMargin));
            xSlabStop = min(xSize, round(i * xSize / splitting + splittingMargin));
            convolvedSlab = uint16(imgaussian(double(stack(xSlabStart:xSlabStop, :, :)), r, r));
            if i == 1
                gaussStack(1:(xSlabStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
            elseif i == splitting
                gaussStack((xSlabStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
            else % i > 1 && i < splitting
                gaussStack((xSlabStart + splittingMargin):(xSlabStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
            end;
            clear convolvedSlab;
        end;
        
        filteredStack = (stack - gaussStack);
        clear gaussStack;
end;
