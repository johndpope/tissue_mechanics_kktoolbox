function mosaic = assign_sinks_to_voxels_dist(L,xs,ys,zs, begin, stop)
%% L is the label matrix of sinks
%% xs, ys and zs contain the sink coordinates that a voxel is associated
%% with.

%%%% the code below can be distributed 
% speed optimization is necessary here

verbose = 1;
if verbose, disp('Assigning sinks to voxels');end

num = max(L(:));    % determines the number of sinks
mosaic = zeros(size(L));
tm = zeros(size(L));
for six = begin:stop,        % loop over the sinks and record the relevant coordinates
    if verbose,disp(['Processing object number: ' num2str(six) ' of ' num2str(num)]);end
    vec = find(L==six);     % gets the indices of the coordinates covered by one sink
    [yv xv zv] = ind2sub(size(L),vec);  % determine the coordinates of the sink valley
    for ix = 1:length(vec)      % loop over the set of coordinates that make up the sink valley
        tm(:) = 0;
        yix = find(ys==yv(ix));
        zix = find(zs==zv(ix));
        tm(xs==xv(ix)) = 1;
        tm(yix) = tm(yix)+1;
        tm(zix) = tm(zix)+1;
        mosaic(tm==3) = six;
    end
end