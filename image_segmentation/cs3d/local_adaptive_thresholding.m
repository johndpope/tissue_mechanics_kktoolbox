function [segmented mosaic] = local_adaptive_thresholding(I, mosaic,num, lbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STEP III : Local adaptive thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1;
use_sct = 0;
dbg = 1;
max_erode = 200;        % maximum while loop iterations for a pair of objects
segmented = zeros(size(I), 'double'); % this will be a label matrix with the objects
temp = zeros(size(I), 'double');
if verbose,disp('Pruning small mosaic regions and performing local adaptive thresholding ...');end
vol = zeros(1,num);  % will store the voxel count for each mosaic region
counter = 0;        % reserved counter for sequential labeling of accepted and segmented objects.
label_index = [];   % saves the index of the objects that are included into the label matrix

%%%%% determine volume lowerbound using kmeans if lbnd not given
if nargin<4||lbnd == 0,
    if verbose,disp(' ');disp('Automatic determination of volume lower bounds using kmeans...');drawnow;end
    for ix = 1:num,
        vol(ix) = sum(mosaic(:)==ix);
    end
    idx = kmeans(vol',2);
    v_large = min(vol(idx==1));v_small = max(vol(idx==2));
    if v_large(1)>v_small(1);v_large = 1; v_small = 2;else v_large = 2;v_small = 1;end %#ok<NASGU>

    volidx = find(vol(v_large));    % these are indices into vol that specify the large volumes
    [vol_sorted] = sort(vol(volidx)); %#ok<FNDSB> % sort these large volumes
    lbnd = vol_sorted(1);
    if verbose,disp(['Lower bound estimate: ' num2str(lbnd)]);disp('Done!');end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ix = 1:num,
    vol(ix) = sum(mosaic(:)==ix);
    if dbg,disp(['object : ' num2str(ix)]);end
    if vol(ix)<lbnd, mosaic(mosaic==ix) = 0;
        if dbg, disp(['eliminated due to small region volume: ' num2str(vol(ix))]);end
    else
        
        [X Y Z] = ind2sub(size(mosaic),find(mosaic==ix));
        block = I(min(X):max(X), min(Y):max(Y), min(Z):max(Z));
        block_mask = mosaic(min(X):max(X), min(Y):max(Y), min(Z):max(Z));
        block(block_mask~=ix) = 0;
        block = mat2gray(block);
        
        if use_sct,
        %%% to do sct
        bb = block(block>0);
        [tt level] = graythresh_sct(bb, 10, 5000);  % use stable count threshold (Russell et al.BPJ 2009)
        block(block<level) = 0;block(block>0) = 1;
        else
        %% to do otsu
        bb = block(block>0);
        level = graythresh(bb); 
        block(block<level) = 0; block(block>0) = 1;
        % uses Otsu's method on the non-zero voxel values only
        end
        
        
        % test whether including this object into the label matrix will
        % create overlap or not. If yes, report and erode till no overlap
        overlap = 0;
        sb = segmented(min(X):max(X), min(Y):max(Y), min(Z):max(Z)); % the region of "segmented" that the block will be added to
        sbvals = unique(sb(:)); % let's see what values we have. Remember these values correspond to  'counter'
        sbvals(sbvals<=0) = [];
        if isempty(sbvals)==0,        % then we should check for potential overlap
            if dbg, disp(['Checking for potential overlap. Values: ' num2str(sbvals(:)')]);end
            for segix = 1:length(sbvals),       % loop over values of sbvals
               val = label_index(sbvals(segix));     %  sbvals correspond to counter values in 'segmented' not object numbers in mosaic
                                                     %  interpret 'counter' values in terms of object number.
               added = sb + double(block*99999);      % the region with both sb and block together
               if any(added(:)==sbvals(segix)+99999) % then we have overlap with the object number "label_index(val)"
                                        % and must erode both objects till we don't have overlap
                  warning('Overlap detected');
                  overlap = 1; %#ok<NASGU>
                  %if verbose, disp('Overlap detected. Eroding both objects.');end
                  segmented(segmented==sbvals(segix)) = 0;    % delete the object from 'segmented' wrt val
                  %%% let's get the old block
                  [Xo Yo Zo] = ind2sub(size(mosaic),find(mosaic==val));
                  blocko = I(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo));
                  block_masko = mosaic(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo));
                  blocko(block_masko~=val) = 0;

                  %%% to do sct
                  if use_sct
                      bbo = blocko(blocko>0);
                      [tt levelo] = graythresh_sct(bbo, 10, 5000);  % use stable count threshold (Russell et al.BPJ 2009)
                      blocko(blocko<levelo) = 0;blocko(blocko>0) = 1;
                  else
                      % to do otsu
                      bbo = blocko(blocko>0);
                      levelo = graythresh(bbo);blocko(blocko<levelo) = 0;blocko(blocko>0) = 1;      % uses Otsu's method on the non-zero voxel values only
                  end
                  
                  wcount = 0;
                  while (any(added(:) == sbvals(segix)+99999) ||wcount>max_erode) % erode both block and blocko till no overlap is detected
                      wcount = wcount + 1;
                      if dbg,disp('Eroding old and new blocks');end
                      seo = strel('ball',3,3);
                      blocko = imerode(blocko,seo);blocko = mat2gray(blocko);     % erode the old shape
                      blocko(blocko>0.5) = 1;blocko(blocko<0.5) = 0;
                      temp(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo)) = ...
                          temp(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo)) + blocko*sbvals(segix);
                      block_old = temp(min(X):max(X), min(Y):max(Y), min(Z):max(Z)); % region of segmented (here temp) that only has the new object
                      
                      block = imerode(block,se); block = mat2gray(block);    % erode the new shape
                      block(block>0) = 1;
                      temp(:) = 0;
                      added = block_old+block*99999;
                  end
                  overlap = 0;
                  %%% add the eroded blocko to 'segmented' ('block' will
                  %%% be added as usual)
                  segmented(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo)) = ...
                    segmented(min(Xo):max(Xo), min(Yo):max(Yo), min(Zo):max(Zo)) + blocko*sbvals(segix);
               end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if overlap==0,
            counter = counter+1;    % increment label counter
            label_index = [label_index ix]; %#ok<AGROW>
            if max(block(:))>1, error('Block values must not exceed 1');end
            segmented(min(X):max(X), min(Y):max(Y), min(Z):max(Z)) = ...
                segmented(min(X):max(X), min(Y):max(Y), min(Z):max(Z)) + double(block*counter);
            if max(segmented(:))>ix, error('Segmented value exceeded current object index.');end
            %if length(label_index)~=max(segmented(:)), error('Discrepancy between label_index and segmented');end
        end
    end
end
if verbose,
    disp(['Number of objects found by local adaptive thresholding: ' num2str(max(segmented(:)))]);
    disp('Done!');
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