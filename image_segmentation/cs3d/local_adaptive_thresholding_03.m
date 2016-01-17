function [segmented mosaic] = local_adaptive_thresholding_03(I, mosaic,num, lbnd, ubnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STEP IV : Local adaptive thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;
dbg = 0;
if nargin<5, ubnd = inf;end
segmented = zeros(size(I), 'double'); % this will be a label matrix with the objects
temp = zeros(size(I), 'double');
if verbose,disp('Performing local adaptive thresholding (light version)...');end
vol = zeros(1,num);  % will store the voxel count for each mosaic region
counter = 0;        % reserved counter for sequential labeling of accepted and segmented objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ix = 1:num,
    
    vol(ix) = sum(mosaic(:)==ix);
    if dbg,disp(['object : ' num2str(ix)]);end
    if vol(ix)<lbnd || vol(ix)>ubnd, mosaic(mosaic==ix) = 0;
        if dbg, disp(['eliminated due to small region volume: ' num2str(vol(ix))]);end
    else
        counter = counter + 1;
%         [X Y Z] = ind2sub(size(mosaic),find(mosaic==ix));
        indx = find(mosaic==ix);
        bb = mat2gray(I(indx));
        %%%  do otsu
        level = graythresh(bb);            % uses Otsu's method on the non-zero voxel values only
        %level = 0.6;
        if verbose, disp(['---------------- threshold level: ' num2str(level)]);end
        sindx = indx(bb>level);
        segmented(sindx) = counter;%segmented(sindx) + counter;
    end
end

%%
% % if verbose,
% %     disp(['Number of objects found by local adaptive thresholding: ' num2str(max(segmented(:)))]);
% %     disp('Done!');
% %     if~isdeployed
% %         h = figure('visible','on');kk_montage(mat2gray(segmented));title('Segmented');
% %         %print -dtiff -r600 intermediate_result_segmented_mosaic.tif;
% %     end
% % end
