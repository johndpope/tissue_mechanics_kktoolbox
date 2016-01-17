function write_j2k_stack(img, filename)
%  write a 3D image file in jpeg 2000 format
% % mkdir('temp');cd temp;
% % if ndims(img) ~= 3
% %     help wr_img
% %     error('Requires image in 3D')
% % end
% % 
% % if nargin < 2
% %     error 'requires image and file name';
% % end
% % if nargin <3, method = 'none';end
% % if strcmpi(method,'jpeg')&&nargin<4, n = 1;end
% % wr_mode = 'overwrite';
% % 
% % 
% % for i=1:size(img,3)
% %     str = sprintf('%s_00%d',filename,i);
% %     warning off;imwrite(img(:,:,i),str, 'jp2', 'CompressionRatio',2.0);warning on
% %     wr_mode = 'append';
% % end
