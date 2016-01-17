function [X_o] = shp_scale(X, scale)
%% uniformly scales the shape X using the vector "scale" with
%% translation
X_o = [];
scale = scale(:);
for ix = 1:size(X,1);
    Xtmp = X(ix,:);
    [xc1 yc1 zc1] = get_xyz_clks(Xtmp);
    
    xc1= xc1.*scale(1);
    yc1= yc1.*scale(2);
    zc1 = zc1.*scale(3);
    
    X_o = [X_o; [xc1(:)' yc1(:)' zc1(:)']];
end

