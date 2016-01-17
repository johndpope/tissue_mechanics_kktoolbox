function indx = closest_point(X,p)
%%% find the index of the point in X that is closes to the coordinates
%%% given in p

dmin = inf;
for ix = 1:length(X)
    d = sqrt(sum((X(ix,:)-p).^2));
    if d<dmin, 
        dmin = d;
        ixmin = ix;
    end
end
indx = ixmin;
