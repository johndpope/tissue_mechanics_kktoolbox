function kk_montage(I,J)
%%%
% if exist('I_display'), I_display(:,:,1,:) = 0;
% else
%     I_display = zeros(size(I,1),size(I,2),1,size(I,4));
% end
warning off;
if nargin ==1
    I_display(:,:,1,:) = I;
%    montage(mat2gray(I_display(:,:,1,:)));
     montage((I_display(:,:,1,:)));
end
if nargin ==2
    I_display(:,:,1,:) = I;
    I_display(:,:,2,:) = 0;
    I_display(:,:,3,:) = J;
    montage((I_display(:,:,:,:)));
end
warning on;