function T = fuse(SHP)
%%% fuse the spherical harmonics objects identified in SHP cell array.
%%% Length SHP is the number of views. angles is a matrix whose rows are
%%% the Euler angles of the view from the microscopy relative to 0 degrees
%%% Author: Khaled Khairy
%%% Uses: rotate_shps, kk_cluster, r_inv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first we need to rotate every "view" back by its angle
% for ix  = 1:length(SHP)
%     ang = angles(ix,:);
%     S = rotate_shps(SHP{ix},-ang, [0 0 0]);
%     SHP{ix} = S;
% end
%% calculate the rotation invariant version of them
disp('Generating the rotation invariance....');
% matlabpool local 8
L_max = 10;
SHP_ri = {};
SHP_ri_p_only = {};
SHP_ri_p_and_o_only = {};
SHP_ang = {};
warning off
parfor vix = 1:length(SHP), % parallel "for" loop over the views
    all_X_o_ri = [];
    all_X_o_ri_p_only = [];
    all_X_o_ri_p_and_o_only = [];
    all_ang = [];
    %disp(['--------------- Processing view: ' num2str(vix) ' of ' num2str(length(SHP))]);
    all_X_o = (SHP{vix});
    for ix = 1:size(all_X_o, 1),    % loop over shapes in this view and get rotation invariance
        %disp(['Processing object [view obj]: ' num2str([vix ix]) ' of ' num2str(size(all_X_o,1))]);
        X_o = all_X_o(ix,1:end);
        X_o = tr(X_o, L_max);
        [X_o X_1 X_2 ang res_p res_o] = r_inv(X_o);
        all_X_o_ri(ix,:) = X_o;
        all_X_o_ri_p_only(ix,:) = X_1;
        all_X_o_ri_p_and_o_only(ix,:) = X_2;
        all_ang(ix,:) = ang;
    end
    SHP_ri{vix} = all_X_o_ri;
    SHP_ri_p_only{vix} = all_X_o_ri_p_only;
    SHP_ri_p_and_o_only{vix} = all_X_o_ri_p_and_o_only;
    SHP_ang{vix} = all_ang;
end
warning on
% matlabpool close
save SHP_ri SHP_ri SHP_ri_p_only SHP_ri_p_and_o_only SHP_ang
disp('Done !');



%% perform shape correspondence
clc;load SHP_ri;    %plot_shp_clouds(SHP_ri);
n = 3;          % the number of neighbors to be considered
[C M Z] = kk_cluster(SHP_ri, n);
save data_Clustering_results C M Z


%% now attach scores to the bins
nbins = max(C);
T_all    = cell2mat(SHP');    % matrix with all shapes view-rotation-matched
T_all_ri = cell2mat(SHP_ri');% matrix with all ri shapes
Tpo = cell2mat(SHP_ri_p_and_o_only');% matrix with all ri shapes
All_ang = cell2mat(SHP_ang');% matrix with all obj angle transformations that would take the shape into its original orientation (not the view orientation)
SC    = cell(nbins,1);
SC_ri = cell(nbins,1);
SA    = cell(nbins,1);
D = zeros(nbins,1);
nObjs = zeros(nbins,1);
[cc IX] = sort(C);  % sorts C and stores the positions in IX
for ix = 1:nbins,  % loop over the bins and calculate score for each
    indx_vec    = IX(cc==ix);
    nobjs(ix)   = length(indx_vec);
    SC{ix}      = T_all(indx_vec,:);
    SC_ri{ix}   = T_all_ri(indx_vec,:);
    SCpo{ix}   = Tpo(indx_vec,:);
    SA{ix}      = All_ang(indx_vec,:);
    D(ix) = 0;
    %%% we need to add all dissimilarities for each bin and average to set a score of
    %%% confidence in the collection within the bin
    for bix = 1:nobjs(ix)
        for bjx = bix+1:nobjs(ix)
            D(ix) = D(ix) + M(indx_vec(bix),indx_vec(bjx));
        end
    end
    D(ix) = D(ix)/nobjs(ix);        % to get the average
end
save
%% FUSE
%%% We have
%%% SC: view-rotation-aligned objects collected in bins, 
%%% SCpo: view-rotation-aligned and parameterization and object-space aligned collected in bins
%%% SA  : Angles that bring the objects from object-aligned back to original rotation around their own axis
%%% cc  : correspondence 
%%% D   : bin quality scores
load
T_av = [];
T = [];
for ix = 1:nbins,
    
%     tt = SC_ri{ix};
%     v1 = tt(1,:);v2 = tt(2,:);v3 = tt(3,:);
    
    vv = sum(SC_ri{ix})/size(SC_ri{ix},1);
    T_av = [T_av;vv];
    ang = SA{ix};   % take the first set of angles
    T(ix,:) = restore_configuration_from_r_inv(T_av(ix,:),ang(3,:));
end
save























