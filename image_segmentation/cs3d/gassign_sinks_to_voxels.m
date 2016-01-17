function L = gassign_sinks_to_voxels(L,xs,ys,zs)
% L is the label matrix of sinks. xs, ys and zs contain the sink 
% coordinates that a voxel is associated with. 
% This means that if we want to know where voxel [a b c] would end
% up, i.e. in which sink, then we need to check [xs(a) ys(b) zs(c)].
% Output: L is preplaced by the mosaic image that reflects the extent of
% influence that a sink has.
% Author: Khaled Khairy, JFRC 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;
GPU = 0;
if verbose && GPU, disp(['GAssigning sinks to voxels (GAS2V).']);end
numSinks = max(L(:));    % determines the number of sinks
v = zeros(numSinks , 3);
sIx = 1:numSinks ;
for six = 1:numSinks ,
    vec = find(L==six);     % gets the indices of the coordinates covered by one sink
    [yv xv zv] = ind2sub(size(L),vec);  % determine the coordinates of the sink valley
    v(six, :) = round([sum(yv) sum(xv) sum(zv)]/length(vec));   % center of mass of each sink
end

if GPU
    %%%%%%%%%%% the process is data-parallel so use GPU
    L = gsingle(L);
    xs = gsingle(xs);
    ys = gsingle(ys);
    zs = gsingle(zs);
end
if verbose && GPU,disp(['GAS2V: Peforming GPU calculation']);end
thresh = 5;

for six = 1:numSinks
    d = sqrt((xs-v(six,2)).^2 + (ys-v(six,1)).^2 + (zs-v(six,3)).^2);
    L(d<thresh) = six;
end

% % gfor ix = 1:numel(L)
% %     L(ix) = 0;
% %     for six = 1:num
% %         dl = sqrt((xs(ix)-v(six,2)).^2 + (ys(ix)-v(six,1)).^2 + (zs(ix)-v(six,3)).^2);
% %         if dl<thresh,L(ix) = six;end
% %     end
% % gend
L = double(L);
