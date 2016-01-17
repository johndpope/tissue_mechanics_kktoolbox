function    [segmented num] = discard_small_volumes(segmented)
% % % discard and classify volumes by classification using kmeans
segmented = bwlabeln(segmented, 6);
disp('Discarding small volumes');
stats = regionprops(segmented,'area');
volume = [];
num = length(stats);
% get the list of volumes
for ix = 1:num,
    volume(ix) = stats(ix).Area;
end
volume = volume(:);


if num>200,
    figure;warning off;
    for ix = 2:2,       %loop over the number of volume categtories if desired (at least 2 for descarding small volumes)
        str = sprintf('idx%d = kmeans(double(volume), %d, ''distance'', ''cityblock'',''emptyaction'', ''singleton'',''replicates'',5);%[silh%d,h] = silhouette(double(volume),idx%d);s(%d) =mean(silh%d); disp([%d mean(silh%d)]);',...
            ix,ix, ix,ix, ix, ix, ix, ix);
        eval(str);
    end
    warning on;
    close;

    % % % discard small volumes based on idx2;
    v_large = min(volume(idx2==1));v_small = max(volume(idx2==2));
    if v_large(1)>v_small(1);v_large = 1; v_small = 2;else v_large = 2;v_small = 1;end
    counter = 0;
    for ix = 1:num,
        if idx2(ix) == v_small,  % if we're looking at a small volume
            segmented(segmented==ix) = 0;
        else
            counter = counter +1;
            segmented(segmented==ix) = counter;
        end
    end
    segmented = bwlabeln(segmented, 6); % just to relabel the remainder
    num = max(segmented(:));
    disp(['After discarding small volumes, the # of objects is: ' num2str(num)]);
    disp('Done!');
    kk_montage(mat2gray(segmented));
end