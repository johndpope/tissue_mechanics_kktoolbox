function mosaic = assign_sinks_to_voxels(L,xs,ys,zs)
% L is the label matrix of sinks. xs, ys and zs contain the sink
% coordinates that a voxel is associated with.
% This means that if we want to know where voxel [a b c] would end
% up, i.e. in which sink, then we need to check [xs(a) ys(b) zs(c)].
%
% Author: Khaled Khairy, JFRC 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verbose = 1;
if verbose, disp(['Assigning sinks to voxels (AS2V).']);end
num = max(L(:));    % determines the number of sinks
mosaic = zeros(size(L));
%tm = cell(num,1);
tm  = zeros(size(L));
for six = 1:num,        % loop over the sinks and record the relevant coordinates
    %tm{six} = zeros(size(L));
    if verbose,disp(['AS2V:Processing sink: ' num2str(six) ' of ' num2str(num)]);end
    vec = find(L==six);     % gets the indices of the coordinates covered by one sink
    [yv xv zv] = ind2sub(size(L),vec);  % determine the coordinates of the sink valley
    if verbose,disp(['AS2V:Found: ' num2str(length(vec)) ' voxels for this sink']);end
    
    for ix = 1:length(vec)      % loop over the set of coordinates that make up the sink valley
        tm(:) = 0;
        yix = find(ys==yv(ix));
        zix = find(zs==zv(ix));
        tm(xs==xv(ix)) = 1;
        tm(yix) = tm(yix)+1;
        tm(zix) = tm(zix)+1;
        mosaic(tm==3) = six;
    end
    
    %     for ix = 1:length(vec)      %% I found this hard to vectorize --- any suggestions?
    %         tm{six}(xs==xv(ix))=1;
    %         tm{six}(ys==yv(ix))=tm{six}(ys==yv(ix))+1;
    %         tm{six}(zs==zv(ix))=tm{six}(zs==zv(ix))+1;
    %     end
end
%
% for six = 1:num,
%     mosaic(tm{six}==3) = six;
% end

