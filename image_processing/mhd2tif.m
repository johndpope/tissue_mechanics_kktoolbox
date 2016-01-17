function D = mhd2tif(filename, flag)
% convert mhd (output from drop3d) to tif file
if nargin==1, flag = 1;end
P = read_mhd(filename);
D = read_raw_data(P);
D = mat2gray(permute(D, [2 1 3]));
if flag,
    write_tif_stack(D, [filename '.tif']);
end;