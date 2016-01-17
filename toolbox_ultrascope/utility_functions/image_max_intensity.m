function Bmx = image_max_intensity(B, dim, plot_flag)
% returns (and plots) the maximum intensity projection of B
if nargin<3, plot_flag = 0;end
Bmx = squeeze(max(B,[],dim));

if plot_flag
imtool(mat2gray(Bmx)); % just to make it look like in Fiji
end