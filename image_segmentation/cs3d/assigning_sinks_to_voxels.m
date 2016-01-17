function mosaic = assigning_sinks_to_voxels(L,num,xs,ys,zs)
%%%% the code below can be distributed 
% speed optimization is necessary here
verbose = 1;
if verbose, disp('Assigning sinks to voxels');end


mosaic = zeros(size(L));
tm = zeros(size(L));
for six = 1:num,        % loop over the sinks and record the relevant coordinates
    if verbose,disp(['Processing object number: ' num2str(six) ' of ' num2str(num)]);end
    vec = find(L==six);
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