function j2c_to_tif()
%% convers all files ending with the extension .j2c to tif files in the current directory
%% deletes the .j2c files as it goes
fns = dir('*.j2c');
n = length(fns);
for ix = 1:n
    fn = fns(ix).name;
    write_stack(imread(fn),[fn '.tif']);
    delete(fn);
end

%%%%%%%%%%%%%%%%%%
function write_stack(img, filename, method,n)
%  write a 3D image file in tif format

% if ndims(img) ~= 3
%     error('Requires image in 3D')
% end

if nargin < 2
    error 'requires image and file name';
end
if nargin <3, method = 'none';end
if strcmpi(method,'jpeg')&&nargin<4, n = 1;end
wr_mode = 'overwrite';
for i=1:size(img,3)
    if strcmpi(method,'jpeg'),
        imwrite(img(:,:,i),filename, 'tif', 'Compression',method,'RowsPerStrip',8*n,'WriteMode', wr_mode);
    else
        imwrite((img(:,:,i)),filename, 'Compression',method,'WriteMode', wr_mode);
    end
    wr_mode = 'append';
end


