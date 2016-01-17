function post_process_cluster_jobs( cluster_path,filename, sub_filenames)

% Creates directories locally if necessary and copies files from the
% cluster to the directories

mkdir shp
mkdir rendered
mkdir segmented
mkdir entropy


%%% copy the SHP files
str = 'shp';
if exist(str)==7, mkdir shp; pause(1);end;
for ix = 1:length(sub_filenames)
    dir_str = sprintf('%s\\Dir_%s',cluster_path,sub_filenames{ix});
    %%% copy the sub image into this directory
    dest = (str);
    file_str = sprintf('data_shp_%s_%s.tif.mat', filename,get_suffix(3,ix));
    src = sprintf('%s\\%s', dir_str, file_str);
    if exist(src)==2,     
        disp(['Copying from: ' src ' to ' dest]);
        copyfile(src,dest,'f');
    end
end


%%% copy the rendered files
str = 'rendered';
if exist(str)==7, mkdir rendered;end;
for ix = 1:length(sub_filenames)
    dir_str = sprintf('%s\\Dir_%s',cluster_path,sub_filenames{ix});
    %%% copy the sub image into this directory
    dest = (str);
    file_str = sprintf('data_render_%s_%s.tif.mat', filename,get_suffix(3,ix));
    src = sprintf('%s\\%s', dir_str, file_str);
    if exist(src)==2,     
        disp(['Copying from: ' src ' to ' dest]);
        copyfile(src,dest,'f');
    end
end

%%% copy the segmentation label files
str = 'segmented';
if exist(str)==7, mkdir segmented;end;
for ix = 1:length(sub_filenames)
    dir_str = sprintf('%s\\Dir_%s',cluster_path,sub_filenames{ix});
    %%% copy the sub image into this directory
    dest = (str);
    file_str = sprintf('data_lat_%s_%s.tif.mat', filename,get_suffix(3,ix));
    src = sprintf('%s\\%s', dir_str, file_str);
    if exist(src)==2,     
        disp(['Copying from: ' src ' to ' dest]);
        copyfile(src,dest,'f');
    end
end

%%% copy the entropy label files
str = 'entropy';
if exist(str)==7, mkdir segmented;end;
for ix = 1:length(sub_filenames)
    dir_str = sprintf('%s\\Dir_%s',cluster_path,sub_filenames{ix});
    %%% copy the sub image into this directory
    dest = (str);
    file_str = sprintf('data_score_%s_%s.tif.mat', filename,get_suffix(3,ix));
    src = sprintf('%s\\%s', dir_str, file_str);
    if exist(src)==2,     
        disp(['Copying from: ' src ' to ' dest]);
        copyfile(src,dest,'f');
    end
end

