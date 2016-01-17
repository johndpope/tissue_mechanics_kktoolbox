function [T E delix ILoc] = score_T(I,T)
% calculate the normalized entropy of the block of image that each object is located
% at and assign that entropy as a contrasting score to the object after
% normalization
E = [];
ILoc = [];
delix = [];
parfor ix = 1:size(T,1),      % loop over the shapes
    disp(['Calculating entropy score for object : ' num2str(ix)]);
    [X]=shp_get_coord(T(ix,1:end-1), 20);
    maxx = round(max(X(:,2))); minx = round(min(X(:,2)));
    maxy = round(max(X(:,1))); miny = round(min(X(:,1)));
    maxz = round(max(X(:,3))); minz = round(min(X(:,3)));
    %%% only score those shapes that lie within the boundary and delete the
    %%% rest
    if minx<1 || maxx>size(I,2) || miny<1 || maxy>size(I,1) || minz<1 ||maxz>size(I,3),
        delix(ix) = 1;
        E(ix) = nan;
        ILoc(ix) = nan
        disp('Object outside boundary ---> deleted');
    else
        delix(ix) = 0;
        Ib = mat2gray(I(miny:maxy, minx:maxx,minz:maxz));
        E(ix) = entropy(Ib);
        ILoc(ix) = mean(Ib(:));
    end
end
T(find(delix==1),:) = [];
E(find(delix==1))   = [];
ILoc(find(delix==1)) = [];
E = E./max(E);  % normalization of Entropy
ILoc = ILoc./max(ILoc);  % normalization of local intensity
