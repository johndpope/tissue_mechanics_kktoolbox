%% this script builds a POVray movie out of a sequence of Ultrascope images
%%% uses : POVray and Fiji
%%% Author: Khaled Khairy
%% configure
clc;
data_root = 'D:\Khaled\Dros_11_01_12_fused_all_movie\';
file_struct = dir([data_root '*.df3']);
source_fn = {};
target_fn = {};
for ix = 1:numel(file_struct),
    source_fn{ix} = file_struct(ix).name;
    target_fn{ix} = [data_root file_struct(ix).name '.pov'];
end


fid = fopen('txt1.txt');
tline = fgetl(fid);
prestr = tline;
while ischar(tline)
    tline = fgetl(fid);
    prestr = [prestr '\n' tline];
end
fclose(fid);

fid = fopen('txt2.txt');
tline = fgetl(fid);
poststr = tline;
while ischar(tline)
    tline = fgetl(fid);
    poststr = [poststr '\n' tline];
end
fclose(fid);

%% generate the pov macro
for ix = 1:numel(file_struct)
    disp(['Generating povray macro : ' target_fn{ix}]);
    [fid msg] = fopen(target_fn{ix}, 'w');
    fprintf(fid, prestr);
    fprintf(fid,'density_file df3 "%s" interpolate 1\n', [data_root source_fn{ix}]);
    fprintf(fid, poststr);
    fclose(fid);
end
