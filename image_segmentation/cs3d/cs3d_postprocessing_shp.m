function [all_T] = cs3d_postprocessing_shp( cluster_path,...
    filename, sub_filenames, margin, wb_thresh, V_max_obj, V_min_obj,...
    A_max_obj, A_min_obj)
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
all_T = [];
for ix = 1:nfiles
    counter = counter + 1;
    
    data_file_str = sprintf('%s\\Dir_%s_%s.tif\\data_shp_%s_%s.tif.mat',...
        cluster_path,filename, get_suffix(3,ix),filename, get_suffix(3,ix));
    disp(['Loading rendered result from fileserver: ' data_file_str]);
    load(data_file_str,'T','T_cleaned', 'wb_vec', 'A_vec', 'V_vec');
    
% %     %% for old version of shp_gen
% %     T_cleaned  = T(find(wb_vec<wb_thresh),:); %#ok<FNDSB>
% %     wb_vec = wb_vec(wb_vec<wb_thresh);
% %     V_vec = V_vec(wb_vec<wb_thresh);
% %     A_vec = A_vec(wb_vec<wb_thresh);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    T = T_cleaned(:,1:end-1);
    
    position_data_file_str = sprintf('data_%s_%s.mat',filename, get_suffix(3,ix));
    disp(['Loading position from local disk: ' position_data_file_str]);
    load(position_data_file_str);
    
    discard = [];
    for oix = 1:size(T,1),  % loop over the objects
         if V_vec(oix)>V_max_obj || A_vec(oix)>A_max_obj ||V_vec(oix)<V_min_obj||A_vec(oix)<A_min_obj,
                          discard = [discard;oix];
         else
        X_o = T(oix,:);
        [xc yc zc] = get_xyz_clks(X_o);
        cm_vec = [xc(1) yc(1) zc(1)]/sqrt(4*pi);
        if (cm_vec(1)<margin || cm_vec(1)> (max(imx(:))- margin) ...
                || cm_vec(2)<margin || cm_vec(2)> (max(imy(:)-margin))), % is the object cm in the overlap region
            discard = [discard;oix];
        else
            % do the translation here
            xc(1) = xc(1)+ (min(width_range(:))-margin)*sqrt(4*pi); %#ok<NODEF>
            yc(1) = yc(1)+ (min(height_range(:))-margin)*sqrt(4*pi); %#ok<NODEF>
            X_o = [xc(:)' yc(:)' zc(:)'];
            T(oix,:) = X_o;
        end
         end
    end
    T(discard(:)',:) =[];
    all_T = [all_T;T];
    disp(size(all_T,1));
end


%% plotting SHP
str_shp = sprintf('save intermediate_result_T_%s.mat T',filename);
h = figure;
fac = 0.8;
plot_shps(shps_scale_individual(all_T,[1 1 1]*fac),20,'random');colorbar;daspect([1 1 1]);
title('Random color coding.');
xlabel('x');ylabel('y');
view(2);
camlight;
saveas(h,[str_shp '.fig']);
print -dtiff -r600 intermediate_results_SHP_energy.tif;

