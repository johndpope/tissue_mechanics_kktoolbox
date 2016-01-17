function r = j2c(fn,flag, compression_ratio)
% minimal code for performing j2c compression and decompression
% "fn"  input file name with complete path
% "flag" settings:  0 ---- read tif file and compress to j2c format
%                   1 ---- read j2c and decompress to tif
%                   10---- read tif file and compress (verbose)
%                   11---- read j2c file and decompress to tif (verbose)
% Author: Khaled Khairy (Keller lab). Copyright HHMI (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% help block
if nargin == 1 && strcmp(fn,'-help'),
   str = sprintf('j2c help:\n Input: >j2c filename flag compression_ratio\nExample: >j2c filename 0 100\n to compress the tif file "filename" 100x\nInput arguments:\nfilename: complete path of file with extension.\n flag:\n 0 === compress, 10 === compress verbose\n 1 === decompress, 11 === decompress verbose\ncompress_ratio: choose compression ratio.\nAuthor: Khaled Khairy (Keller lab) JFRC (HHMI) 2010.');
   disp(str);
   r = 1;
   return
end
r = 0;
if nargin<2, disp('Two input arguments required. Example: j2c("filename.tif", 0).Type j2c -help for more info.');r = -1;return;end
if nargin<3, compression_ratio = 100;end

if isdeployed, 
    flag = str2num(flag);
    if strcmp(class(compression_ratio),'char'),compression_ratio = str2num(compression_ratio);end
end

if flag ==0,    % perform compression, i.e. convert tif->j2c
    [X im_info] = read_stack(fn);
    if compression_ratio==1,
        imwrite(X,[fn '.j2c'],'Mode', 'lossless', 'Comment', im_info);
    else
        imwrite(X,[fn '.j2c'],'CompressionRatio', compression_ratio, 'Comment', im_info);
    end
elseif flag ==1 % perform decompression, i.e. conver j2c->tif
    write_stack(imread(fn),[fn '.tif']);
elseif flag ==10,    % perform compression, i.e. convert tif->j2c
    tic;disp(['Reading from disk ... ' fn]);
    [X im_info] = read_stack(fn);
    toc
    disp(['Image is type: ' class(X)]);
    
    tic;disp('Comressing and writing to disc');
    if compression_ratio==1,
        imwrite(X,[fn '.j2c'],'Mode', 'lossless', 'Comment', im_info);
    else
        imwrite(X,[fn '.j2c'],'CompressionRatio', compression_ratio, 'Comment', im_info);
    end
    toc
    
elseif flag ==11 % perform decompression, i.e. conver j2c->tif
    tic;disp(['Decompressing and reading from disk ... ' fn]);
    X = imread(fn);
    toc
    disp(['Image is type: ' class(X)]);
    
    tic;disp('Writing to disc');
    write_stack(X,[fn '.tif']);
    toc
else
    disp('Could not recognize input flag. Aborting!');r = -1;return;
end
end


function [I im_info] = read_stack(filename)
%%% load a tif stack file
info = imfinfo(filename);%disp(info);
zdim = length(info);
im = imread(filename,1);
im_info = cell(zdim,1);
if size(im,3)==1,
    I = zeros([size(im,1), size(im,2), zdim], class(im));
    for ix = 1:zdim
        fn = fieldnames(info(ix));
        sc = struct2cell(info(ix));
        str = '';
        for jx = 1:length(sc),
            if ~strcmp(class(sc{jx}), 'char'), sc{jx}= num2str(sc{jx});end
            str = [str [fn{jx} ' ' sc{jx} '\n']];
        end
        im_info{ix} = str;
        I(:,:,ix) = imread(filename,ix);
    end
else
    I = zeros([size(im,1), size(im,2), 3], class(im));
    for ix = 1:zdim
        fn = fieldnames(info(ix));
        sc = struct2cell(info(ix));
        im_info{ix} = [fn,sc];
        I(:,:,ix) = imread(filename,ix);
    end
end
end

function write_stack(img, filename, method,n)
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
        imwrite(img(:,:,i),filename, 'tif', 'Compression',method,'RowsPerStrip',8*n,'WriteMode', wr_mode);
    else
        imwrite((img(:,:,i)),filename, 'Compression',method,'WriteMode', wr_mode);
    end
    wr_mode = 'append';
end
end


