
clc;
fiji = '"C:\Khaled_code\Fiji 091014 x64\fiji-win64.exe"';
macro_name = 'prepare_image_stackCM1.ijm';
data_root = 'Y:\\Shared\\Processing\\Ultrascope 1\\11-01-10\\Drosophila_Embryo_Histone-eGFP_10x_2µm_Bidirectional_DualCamera_TL_35s.fused\\';
contrast_slice = 64;
file_struct = dir(['Y:\Shared\Processing\Ultrascope 1\11-01-10\Drosophila_Embryo_Histone-eGFP_10x_2µm_Bidirectional_DualCamera_TL_35s.fused\' '*CM0*.tif']);
source_fn = {};
target_fn = {};
for ix = 1:numel(file_struct),
    source_fn{ix} = file_struct(ix).name;
    target_fn{ix} = ['Y:\\Shared\\Processing\\Ultrascope 1\\11-01-10\\CM0\\CM0_out_' num2str(ix) '.tif'];
end

%% generate the fiji macro
[fid msg] = fopen(macro_name, 'w');
for ix = 1:numel(file_struct)
    fprintf(fid,'open("%s%s");\n', data_root, source_fn{ix});
    %fprintf(fid,'run("Flip Vertically", "stack");');
    fprintf(fid,'makeRectangle(328, 446, 1560, 614);run("Crop");\n');
    fprintf(fid,'setSlice(%d);\n', contrast_slice);
    fprintf(fid,'run("Enhance Contrast", "saturated=0.35");\n');
    %fprintf(fid,'run("Despeckle", "stack");');
    fprintf(fid,'run("Subtract Background...", "rolling=50 stack");\n');
    fprintf(fid,'setSlice(%d);\n', contrast_slice);
    fprintf(fid,'run("Enhance Contrast", "saturated=0.35");\n');
    fprintf(fid,'run("Z Project...", "start=1 stop=135 projection=[Max Intensity]");');
    fprintf(fid,'saveAs("Tiff", "%s");\n',...
        data_root, target_fn{ix});
    fprintf(fid,'run("Close All");\n');
end
fprintf(fid, 'run(''Quit'');');
fclose(fid);

%% call fiji
str = sprintf('%s -macro %s',fiji, [pwd '\' macro_name]);
system(str);