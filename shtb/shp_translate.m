function [X_o] = shp_translate(X, pos)
%% set new translation position according to pos
% X = X(:)';
X_o = [];
pos = pos*sqrt(4*pi);%%% to obtain the translation equivalent for the clks, pos needs to be multiplied by sqrt(4*pi)
for ix = 1:size(X,1);
    Xtmp = X(ix,:);
    [xc2 yc2 zc2] = get_xyz_clks(Xtmp);
    xc2(1)= pos(ix,1);
    yc2(1) = pos(ix,2);
    zc2(1) = pos(ix,3);
    X_o = [X_o; [xc2(:)' yc2(:)' zc2(:)']];
end