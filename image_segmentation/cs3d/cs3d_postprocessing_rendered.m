function [Xc, Fc] = cs3d_postprocessing_rendered(...
    cluster_path,...
    filename,...
    sub_filenames)
% Collect the segmented individual data files and fuse them into one
% dataset.
% What we need is:
%           - The position of the segmented sub image (obtained from local directory)
%           - The path to the segmented rendered data with all_X and all_F
%           for each sub-image
%           - Translational correction on all_X obtained from sub-images.
%           - Gathering translated all_X into an overall all_X for all
%           objects deleting objects in regions of overlap from the next
%           region, i.e. objects from the first sub-image in the sequence
%           have precedence over once from the next sub-image.
%           - Deletion of small objects.
%           - optional cleaning of objects and SH parameterization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfiles = length(sub_filenames);
counter = 0;
Xc = [];
Fc = [];
for ix = 1:nfiles
    counter = counter + 1;
    
    rendered_data_file_str = sprintf('%s\\Dir_%s_%s.tif\\data_render_%s_%s.tif.mat',...
        cluster_path,filename, get_suffix(3,counter),filename, get_suffix(3,counter));
    disp(['Loading rendered result from fileserver: ' rendered_data_file_str]);
    load(rendered_data_file_str,'all_X','all_F');
    
    position_data_file_str = sprintf('data_%s_%s.mat',filename, get_suffix(3,counter));
    disp(['Loading position from local disk: ' position_data_file_str]);
    load(position_data_file_str);
    
    %%% translate the coordinates according to the data stored in
    %%% height_range and width_range
    cm_vec = [];
    for oix = 1:length(all_X),  % loop over the objects
        X = all_X{oix};
        cm_vec = [sum(X,1)./length(X)];
        if (cm_vec(1)<10 || cm_vec(2)<10), % is the object cm in the overlap region
            all_X(oix) = [];
        else
            % do the translation here
            X = [(X(:,1)+ max(width_range(:))-10) (X(:,2)+ max(height_range(:))-10) X(:,3)];
            all_X{oix} = X;
        end
    end
    Xc = [Xc all_X];
    Fc = [Fc all_F];
    
end
%% all objects must go through a rendering cleaning using iso2mesh
%% functions first
global ISO2MESH_SESSION; ISO2MESH_SESSION= filename;
global ISO2MESH_TEMP;ISO2MESH_TEMP = './';
disp('-----------------------------------------------------------');
disp(['---------------- Cleaning ' num2str(length(Xc)) ' surfaces. Please be patient.']);
tic;[Xc, Fc] = clean_all_XF(Xc, Fc, 30, 0);
str = sprintf('save data_render_clean_%s.mat Xc Xc;', filename);eval(str);
toc;disp('Surface triangulation cleaning completed !!');
%% Plotting and saving
h = figure(1); set(h,'visible','off');
plot_all_X(Xc,Fc);view(3);set(h,'visible','on');camlight;drawnow;
str = sprintf('print -dtiff -r600 intermediate_result_cleaned_%s.tif',filename);eval(str);
str  = sprintf('intermediate_result_cleaned_%s.fig',filename);
saveas(h,str);
str = sprintf('save data_intermediate_rendered_cleaned_%s.mat Xc Fc filename',filename);eval(str);

%% SHP representation
tic;
[T, T_cleaned, V_vec, A_vec, wb_vec, k_g_vec, h_vec]  = shp_gen(...
     Xc, Fc, 10, 30, 100);
str_shp = sprintf('data_shp_%s.mat',filename);
str = sprintf('save data_shp_%s.mat T T_cleaned V_vec A_vec wb_vec k_g_vec h_vec;',...
              filename);eval(str);
toc;
disp('SHP fitting completed !!');

%% plotting SHP
load(str_shp);
h = figure('visible','off');
plot_shp_shapes(T_cleaned(:,1:end-1),20,'blue');
title('Random color coding.');
saveas(h,[str_shp '.fig']);
%print -dtiff -r600 intermediate_results_SHP_energy.tif;
clf
close(h);
