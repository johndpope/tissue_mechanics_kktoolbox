function write_tif_stack(img, filename, method,n)
%  write a 3D image file in tif format

if ndims(img) ~= 3
    error('Requires image in 3D')
end

if nargin < 2
    error 'requires image and file name';
end
if nargin <3, method = 'none';end
if strcmpi(method,'jpeg')&&nargin<4, n = 1;end
wr_mode = 'overwrite';
for i=1:size(img,3)
    if strcmpi(method,'jpeg'),
    warning off
    imwrite(img(:,:,i),filename, 'tif', 'Compression',method,'RowsPerStrip',8*n,'WriteMode', wr_mode);
    warning on
    else
    warning off
    imwrite((img(:,:,i)),filename, 'Compression',method,'WriteMode', wr_mode);
   
warning on
    end
    wr_mode = 'append';
end
