function [segmented mosaic maxcoord] = local_adaptive_thresholding_02(I, mosaic,num, lbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STEP III : Local adaptive thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1;
use_sct = 0;
dbg = 1;
max_erode = 200;        % maximum while loop iterations for a pair of objects
segmented = zeros(size(I), 'double'); % this will be a label matrix with the objects
temp = zeros(size(I), 'double');
if verbose,disp('Performing local adaptive thresholding (light version)...');end
vol = zeros(1,num);  % will store the voxel count for each mosaic region
E = zeros(1,num);
counter = 0;        % reserved counter for sequential labeling of accepted and segmented objects.
label_index = [];   % saves the index of the objects that are included into the label matrix
maxcoord = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ix = 1:num,
    counter = counter + 1;
    vol(ix) = sum(mosaic(:)==ix);
    if dbg,disp(['object : ' num2str(ix)]);end
    if vol(ix)<lbnd, mosaic(mosaic==ix) = 0;
        if dbg, disp(['eliminated due to small region volume: ' num2str(vol(ix))]);end
    else
        [X Y Z] = ind2sub(size(mosaic),find(mosaic==ix));
        indx = find(mosaic==ix);
        
        Iblock = I(min(X):max(X), min(Y):max(Y), min(Z):max(Z));
        block = mat2gray(Iblock);
        
        %%% check block entropy
        E(ix) = entropy(block);disp(['---------------- Block entropy:  ' num2str(E(ix))]);
        figure(6);kk_montage(mat2gray(Iblock));
        if use_sct,
            %%%  do sct
            bb = block(block>0);
            [level tt] = graythresh_sct(bb, 10, 5000);  % use stable count threshold (Russell et al.BPJ 2009)
            block(block<level) = 0;block(block>0) = 1;
        else
            %%%  do otsu
            bb = block(block>0);
            level = graythresh(bb);
            disp(['---------------- threshold level: ' num2str(level)]);
            block(block<level) = 0; block(block>0) = 1;
            % uses Otsu's method on the non-zero voxel values only
        end
        figure(8);kk_montage(block);
        mIblock = medfilt3(mat2gray(Iblock),5);
        [d1 d2 d3] = ind2sub(size(block), find(mIblock == max(mIblock(:))));
        maxcoord= [maxcoord; d1(1) d2(1) d3(1)];
        segmented(min(X):max(X), min(Y):max(Y), min(Z):max(Z)) = ...
            segmented(min(X):max(X), min(Y):max(Y), min(Z):max(Z)) + double(block*counter);
        if max(segmented(:))>ix, error('Segmented value exceeded current object index.');end
        %if length(label_index)~=max(segmented(:)), error('Discrepancy between label_index and segmented');end
        figure(7);kk_montage(mat2gray(segmented));
    end
end
if verbose,
    disp(['Number of objects found by local adaptive thresholding: ' num2str(max(segmented(:)))]);
    disp('Done!');
    figure;hist(double(E));
    if~isdeployed
        h = figure('visible','off');kk_montage(mat2gray(segmented));
        print -dtiff -r600 intermediate_result_segmented_mosaic.tif;
    end
end

%%% Note on Otsu's method:
%% The algorithm assumes that the image to be thresholded contains two
%% classes of pixels (e.g. foreground and background) then calculates the
%% optimum threshold separating those two classes so that their combined
%% spread (intra-class variance) is minimal.