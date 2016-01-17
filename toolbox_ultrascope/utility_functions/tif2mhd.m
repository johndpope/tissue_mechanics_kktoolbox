function D = tif2mhd(filename)
% convert tif image to metafile mhd format readable by dropreg3d.exe
info = imfinfo(filename);
xdim = info(1).Width;
ydim = info(2).Height;
zdim = size(info,1);

fid = fopen([filename '.mhd'],'w');
fprintf(fid,'ObjectType = Image\nNDims = 3\nBinaryData = True\nDimSize = %d %d %d\nElementSize = 1.0    1.0     1.0\nElementType = MET_UCHAR\nElementDataFile = %s',...
             xdim, ydim, zdim, [filename '.raw']);
fclose(fid);

P = read_mhd([filename '.mhd']);
D = read_tif_stack(filename);
D = mat2gray(permute(D, [2 1 3]))*255; % convert to 8bit
%D = (permute(D, [2 1 3])); % convert to 8bit

write_mhd(D,P);

         