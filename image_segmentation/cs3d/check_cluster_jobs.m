function finished = check_cluster_jobs(cluster_path, filename, sub_filenames)  
% - checks whether all output has been produced and lists the missing jobs
verbose = 0;
finished = 1;
counter = 0;
for ix = 1:length(sub_filenames)
    sub_fn = sprintf('%s_%s.tif', filename, get_suffix(3,ix));
    file_str = sprintf('data_score_%s_%s.tif.mat', filename,get_suffix(3,ix));
    file_path = sprintf('%s\\Dir_%s\\%s',cluster_path,sub_fn, file_str);
   
    if exist(file_path)==2, 
        if verbose, disp([file_path ' exists']);end;
    else
        if verbose, disp([file_path ' missing']);end
        finished = 0;
        counter = counter + 1;
    end
end
disp([num2str((length(sub_filenames)-counter)/length(sub_filenames)*100) ' % of jobs  finished !']);







