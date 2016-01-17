
clc;
fiji = '"C:\Khaled_code\Fiji 091014 x64\fiji-win64.exe"';
macro_name = 'prepare_4D_movie.ijm';
data_root = 'D:\\Khaled\\Dros_11_01_12_fused_all\\';
contrast_slice = 64;
file_struct = dir(['D:\Khaled\Dros_11_01_12_fused_all\' '*.tif']);
source_fn = {};
target_fn = {};
for ix = 1:numel(file_struct),
    source_fn{ix} = file_struct(ix).name;
    target_fn{ix} = ['D:\\Khaled\\Dros_11_01_12_fused_all_movie_02\\out_' num2str(ix) '.tif'];
end

%% generate the fiji macro
step = 1;
[fid msg] = fopen(macro_name, 'w');
for ix = 1:step:numel(file_struct)
    angle = ix;disp(angle);
    fprintf(fid,'open("%s%s");\n', data_root, source_fn{ix});
    fprintf(fid,'run("Flip Horizontally", "stack");');
    %fprintf(fid,'run("Properties...", "channels=1 slices=116 frames=1 unit=pixel pixel_width=0.363 pixel_height=0.363 voxel_depth=2 frame=[0 sec] origin=0,0");\n');
    fprintf(fid,'setSlice(%d);\n', contrast_slice);
    fprintf(fid, 'resetMinAndMax();run("8-bit");');
    %fprintf(fid,'run("Enhance Contrast", "saturated=0.35");\n');
    %fprintf(fid,'run("Despeckle", "stack");');
    fprintf(fid,'run("TransformJ Scale", "x-factor=1.0 y-factor=1.0 z-factor=5.50964 interpolation=[cubic B-spline]");');
    fprintf(fid,'run("Properties...", "channels=1 slices=639 frames=1 unit=pixel pixel_width=1.0000000 pixel_height=1.0000000 voxel_depth=1 frame=[0 sec] origin=0,0");');
    fprintf(fid,'run("Subtract Background...", "rolling=100 stack");\n');
    fprintf(fid,'setSlice(%d);\n', contrast_slice*5);
    fprintf(fid,'run("Enhance Contrast", "saturated=0.35");\n');

    fprintf(fid, 'run("TransformJ Rotate", "z-angle=0.0 y-angle=0.0 x-angle=%d interpolation=[cubic B-spline] background=0.0");', angle);
    fprintf(fid,'run("DF3 ...", "save=[%s.df3] choose_data_output_format=[8-bit unsigned integer] stretch_contrast");', target_fn{ix});
    fprintf(fid,'run("Z Project...", "start=1 stop=639 projection=[Max Intensity]");');
    fprintf(fid,'saveAs("Tiff", "%s");\n',target_fn{ix});
    fprintf(fid,'run("Close All");\n');
end
fprintf(fid, 'run(''Quit'');');
fclose(fid);

%% call fiji
str = sprintf('%s -macro %s',fiji, [pwd '\' macro_name]);
system(str);