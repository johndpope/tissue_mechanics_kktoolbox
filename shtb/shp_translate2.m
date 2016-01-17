function [X_o] = shp_translate2(X, pos)
%% set new translation position according to pos
%% to obtain the translation equivalent for the clks, pos needs to be multiplied by
%% sqrt(4*pi)
% X = X(:)';
X_o = zeros(size(X));
pos = pos*sqrt(4*pi);
for ix = 1:size(X,1);
    Xtmp = X(ix,:);
    [xc2 yc2 zc2] = get_xyz_clks(Xtmp);
    xc2(1)= xc2(1) + pos(1);
    yc2(1) =yc2(1) +  pos(2);
    zc2(1) =zc2(1) +  pos(3);
    X_o(ix,:) = [[xc2(:)' yc2(:)' zc2(:)']];
end