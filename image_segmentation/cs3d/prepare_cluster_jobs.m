function exit_flag = prepare_cluster_jobs(...
    cluster_path,...                                    % generally cs3d_runs
    cs3d_path,...                                       % generally cs3d_compiled
    run_libraries_path,...                              % generally g/mcmc/khaled (UNIX style)
    run_job_path,...                                    % generally /g/mcmc/khaled/cs3d_runs//the filename (UNIX style)
    matlab_libraries_path,...                           % generally /g/mcmc/khaled/matlab_libraries/v76 (UNIX style)
    compiled_cs3d_command, ...                          % the cs3d executable with path (UNIX style)
    dsfac, zfac,...                                     %downsampling factor and interpolation in z 
    distributed,sub_filenames,...                            % general configuration
    aniso_iter,delta_t,kappa,option, vx, vy, vz,...  % for anisotropic diffusion
    dgvf_iter, miu, dt, dx, dy, dz,...                  % for dgvf configuration
    gft_niter, field_thresh, ...                        % configure gft
    lbnd,...                                            % configure local adaptive thresholding
    V_max, V_min,...                                    % configure rendering
    laplace_iter,L_max, gdim, newton_iter)

% - Creates directories on the project space (fileserver) path given by
%   cluster_path. 
% - Copies into the appropriate directories the image files.
% - Creates a set of scripts for executing cs3d on each image file.


%% start the job_list_script and create the directories
mkdir(cluster_path);        % create general directory for organizing the data processing and results for the upper level
fn_str = sprintf('%s\\cs3d_array.script',cluster_path);
fid_job_list_script = fopen(fn_str, 'w');                  % start a job list script
text_str = sprintf('#!/bin/bash');fwrite(fid_job_list_script, text_str, 'char');

for ix = 1:length(sub_filenames)
    %%% create the directory for organizing processing of the sub images
    dir_str = sprintf('%s\\Dir_%s',cluster_path,sub_filenames{ix});
    disp(['Creating directory: ' dir_str]);
    mkdir(dir_str);
    
    %%% copy the sub image into this directory
    src = sub_filenames{ix};
    dest = sprintf('%s\\%s', dir_str, sub_filenames{ix});
    disp(['Copying from: ' src ' to ' dest]);
    copyfile(src,dest,'f')
    
    %%% write the script functions with parameters and paths for execution
    fn_str = sprintf('%s\\%s.script',cluster_path,sub_filenames{ix});
    fid = fopen(fn_str, 'w');
%     text_str1 = sprintf('#!/bin/bash\n#PBS -l ncpus=1\n#PBS -j oe\n#PBS -m ea\n# Check on some basics:\necho "Running on host: " `hostname`\necho "Changing to directory from which script was submitted."\n#cd $PBS_O_WORKDIR\necho "Current working directory is now: " `pwd`\n# print some environment variables FYI\necho "This job was submitted by user: $PBS_O_LOGNAME"\necho "This job was submitted to host: $PBS_O_HOST"\necho "This job was submitted to queue: $PBS_O_QUEUE"\necho "PBS working directory: $PBS_O_WORKDIR"\necho "PBS job id: $PBS_JOBID"\necho "PBS job name: $PBS_JOBNAME"\necho "PBS environment: $PBS_ENVIRONMENT"\necho " "\necho "This script is running on the PBS MOS node `hostname` "\necho "Job started on: " `date`\n');
    text_str1 = sprintf('#!/bin/bash\n# Check on some basics:\necho "Running on host: " `hostname`\necho "Changing to directory from which script was submitted."\n#cd $PBS_O_WORKDIR\necho "Current working directory is now: " `pwd`\n# print some environment variables FYI\necho "This job was submitted by user: $PBS_O_LOGNAME"\necho "This job was submitted to host: $PBS_O_HOST"\necho "This job was submitted to queue: $PBS_O_QUEUE"\necho "PBS working directory: $PBS_O_WORKDIR"\necho "PBS job id: $PBS_JOBID"\necho "PBS job name: $PBS_JOBNAME"\necho "PBS environment: $PBS_ENVIRONMENT"\necho " "\necho "This script is running on the PBS MOS node `hostname` "\necho "Job started on: " `date`\n');

    run_sub_job_path = sprintf('%sDir_%s',run_job_path, sub_filenames{ix});
    %text_str2 = sprintf('cd %s\n./run_libraries_76.sh\ncd %s',run_libraries_path,run_sub_job_path);
    text_str2 = sprintf('cd %s',run_sub_job_path);
    fwrite(fid, text_str1, 'char');
    fwrite(fid, text_str2, 'char');
   
    text_str3 = sprintf('\nlibrary_path="%s"	# where are the matlab libraries ?',matlab_libraries_path);fwrite(fid, text_str3, 'char');
    text_str = sprintf('\ndsfac="%.4f"					# downsampling factor --- leave 1 for now',dsfac);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nzfac="%.4f"					# z resolution factor (to make isotropic)-- leave 1 for now',zfac);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndistributed="%d"      # set to zero if Distributed toolbox is not installed or when using cluster',distributed);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nfilename="%s"				# image data',sub_filenames{ix});fwrite(fid, text_str, 'char');
    text_str = sprintf('\naniso_iter="%d"             # # of anisotropic diffusion iterations (zero is also ok)',aniso_iter);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndelta_t="%.6f"             # for meaning of remaining anisotropic diffusion parameters please',delta_t);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nkappa="%.2f"                 # type >help anisodiff3D. Also see:',kappa);fwrite(fid, text_str, 'char');
    text_str = sprintf('\noption="%d"                 # Perona and Malik 1990.',option);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nvxsx="%.4f"					# voxel spacing x',vx);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nvxsy="%.4f"					# voxel spacing y',vy);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nvxsz="%.4f"					# voxel spacing z',vz);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndgvf_iter="%d"             # number of diffusion gradient iterations',dgvf_iter);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nmiu="%.2f"                    # large values make vector field smoother (good value: 5)',miu);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndt="%.6f"                 # try to keep dt small (smaller than 0.01 for instance))',dt);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndx="%.2f"                     # don''t change this unless you know what your doing',dx);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndy="%.2f"                     # don''t change this unless you know what your doing',dy);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ndz="%.2f"                     # don''t change this unless you know what your doing',dz);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nfield_thresh="%.4f"       # threshold for points in the vector field to consider to propagate (zero is ok but slower)',field_thresh);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ngft_niter="%d"      # number of voxels estimated between furthest voxel from its sink',gft_niter);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nlbnd="%d"          # lowerbound on the size of mosaic region to process',lbnd);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nV_max="%d"        # Maximum volume of object to consider for surface rendering',V_max);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nV_min="%d"         # Minimum volume of object to consider for surface rendering',V_min);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nlaplace_iter="%d"  # smoothing of the rendered surface triangulations prior to SHP fitting.',laplace_iter);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nL_max="%d"          # maximum L order for SHP fitting. Increase for complicated shapes.',L_max);fwrite(fid, text_str, 'char');
    text_str = sprintf('\ngdim="%d"          # mesh size for L fitting (20->60 is good)',gdim);fwrite(fid, text_str, 'char');
    text_str = sprintf('\nnewton_iter="100"   # number of fitting iterations. Increase for complicated shapes(30->100 is good)',newton_iter);fwrite(fid, text_str, 'char');
    text_str = sprintf('\n%s $library_path $dsfac $zfac $distributed $filename $aniso_iter $delta_t $kappa $option $vxsx $vxsy $vxsz $dgvf_iter $miu $dt $dx $dy $dz $gft_niter $field_thresh $lbnd $V_max $V_min $laplace_iter $L_max $gdim $newton_iter',compiled_cs3d_command);
    fwrite(fid, text_str, 'char');
    fclose(fid);
    
    %%% fill in the line for this sub image
    script_name = sprintf('%s%s.script', run_job_path,sub_filenames{ix} );
    text_str = sprintf('\n/usr/pbs/bin/qsub  -l select=1:ncpus=2,mppmem=4000Mb -q clusterng %s',script_name);fwrite(fid_job_list_script, text_str, 'char');
end
fclose(fid_job_list_script);
exit_flag = 1;





