function [S] = shps_scale_individual(T, scale)
%% uniformly scales the shapes T using scale, without scaling relative distances
%% Input: T is an nxm matrix with n shapes and m is the coefficient vector

if mod(size(T,2),3)==1, T = T(:,1:end-1);end    % just to make sure that we don't have that extra column
S = [];
if length(scale(:))==1, scale = [scale scale scale];end

for ix = 1:size(T,1);
    [xc1 yc1 zc1] = get_xyz_clks(T(ix,:));
    xc = xc1(1:end)*scale(1);xc(1) = xc1(1);
    yc = yc1(1:end)*scale(2);yc(1) = yc1(1);
    zc = zc1(1:end)*scale(3);zc(1) = zc1(1);
 
    S = [S; [xc(:)' yc(:)' zc(:)']];
end
