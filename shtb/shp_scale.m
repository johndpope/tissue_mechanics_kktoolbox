function [S] = shp_scale(T, scale)
%% uniformly scales the shape X using scale 
%% Input: T is an nxm matrix with n shapes and m is the coefficient vector

%if mod(size(T,2),3)==1, T = T(:,1:end-1);end    % just to make sure that we don't have that extra column
S = [];
if length(scale(:))==1, scale = [scale scale scale];end

for ix = 1:size(T,1);
    [xc yc zc] = get_xyz_clks(T(ix,:));
    [xc1 yc1 zc1] = get_xyz_clks(T(ix,:));
    xc(2:end)= xc1(2:end)*scale(1);
    yc(2:end) = yc1(2:end)*scale(2);
    zc(2:end) = zc1(2:end)*scale(3);
    S = [S; [xc(:)' yc(:)' zc(:)']];
end
