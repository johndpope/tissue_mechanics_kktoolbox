function T = fuse_02(SHP,n)
%%% fuse the spherical harmonics objects identified in SHP cell array.
%%% Length SHP is the number of views. angles is a matrix whose rows are
%%% the Euler angles of the view from the microscopy relative to 0 degrees
%%% n is the number of neighbors to be considered
%%% Author: Khaled Khairy
%%% Uses: rotate_shps, kk_cluster, r_inv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2,n = 4;end
T = [];
%% calculate the rotation invariant version of them
disp('Generating the rotation invariance....');
% matlabpool local 8
L_max = get_L_max(SHP{1}(1,:));
all_X_o = cell2mat(SHP');
all_X_o_ri = zeros(size(all_X_o));
all_X_o_ri_p_only = zeros(size(all_X_o));
all_X_o_ri_p_and_o_only = zeros(size(all_X_o));
all_ang = zeros(size(all_X_o, 1), 6);

parfor ix = 1:size(all_X_o,1), %
    %disp(['Processing object [view obj]: ' num2str([vix ix]) ' of ' num2str(size(all_X_o,1))]);
    X_o = all_X_o(ix,1:end);
    X_o = tr(X_o, L_max);
    [X_o X_1 X_2 ang res_p res_o] = r_inv(X_o);
    all_X_o_ri(ix,:) = X_o;
    all_X_o_ri_p_only(ix,:) = X_1;
    all_X_o_ri_p_and_o_only(ix,:) = X_2;
    all_ang(ix,:) = ang;
end

save all_X_o_ri

%%% reorganize the objects into view data structures
SHP_ri = {};
SHP_ri_p_only = {};
SHP_ri_p_and_o_only = {};
SHP_ang = {};
pos = 1;
for vix = 1:length(SHP)
    vec = pos:pos-1+ size(SHP{vix},1)'
    SHP_ri{vix} = all_X_o_ri(vec, :);
    SHP_ri_p_only{vix} = all_X_o_ri_p_only(vec, :);
    SHP_ri_p_and_o_only{vix} = all_X_o_ri_p_and_o_only(vec, :);
    SHP_ang{vix} = all_ang(vec,:);
    pos = pos+ size(SHP{vix},1);
end


% matlabpool close
save SHP_ri SHP_ri SHP_ri_p_only SHP_ri_p_and_o_only SHP_ang
disp('Done !');
%% look at the resulting shape invariances
X1 = SHP_ri_p_only{1};
X2 = SHP_ri_p_only{2};
X3 = SHP_ri_p_only{3};

A1 = SHP_ang{1};
A2 = SHP_ang{2};
A3 = SHP_ang{3};

% % 
% % X11 = X1(1,:); plot_shps(X11);camlight
% % X12 = X1(2,:); plot_shps(X12);camlight
% % 

X13 = X1(3,:); figure;plot_shps(X13);axis equal;camlight
X23 = X2(3,:); figure;plot_shps(X23);axis equal;camlight
X33 = X3(3,:); figure;plot_shps(X33);axis equal;camlight
% % 
% % plot_shp_clouds(SHP_ri');
%% perform shape correspondence
disp('Performing clustering operation');
load SHP_ri; 
n = 3;
[C M Z] = kk_cluster(SHP_ri_p_only, n);

save data_Clustering_results C M Z



%%
% % % % %% now attach scores to the bins
% % % % nbins = max(C);
% % % % T_all    = cell2mat(SHP');    % matrix with all shapes view-rotation-matched
% % % % T_all_ri = cell2mat(SHP_ri');% matrix with all ri shapes
% % % % Tpo = cell2mat(SHP_ri_p_and_o_only');% matrix with all ri shapes
% % % % All_ang = cell2mat(SHP_ang');% matrix with all obj angle transformations that would take the shape into its original orientation (not the view orientation)
% % % % SC    = cell(nbins,1);
% % % % SC_ri = cell(nbins,1);
% % % % SA    = cell(nbins,1);
% % % % D = zeros(nbins,1);
% % % % nObjs = zeros(nbins,1);
% % % % [cc IX] = sort(C);  % sorts C and stores the positions in IX
% % % % for ix = 1:nbins,  % loop over the bins and calculate score for each
% % % %     indx_vec    = IX(cc==ix);
% % % %     nobjs(ix)   = length(indx_vec);
% % % %     SC{ix}      = T_all(indx_vec,:);
% % % %     SC_ri{ix}   = T_all_ri(indx_vec,:);
% % % %     SCpo{ix}   = Tpo(indx_vec,:);
% % % %     SA{ix}      = All_ang(indx_vec,:);
% % % %     D(ix) = 0;
% % % %     %%% we need to add all dissimilarities for each bin and average to set a score of
% % % %     %%% confidence in the collection within the bin
% % % %     for bix = 1:nobjs(ix)
% % % %         for bjx = bix+1:nobjs(ix)
% % % %             D(ix) = D(ix) + M(indx_vec(bix),indx_vec(bjx));
% % % %         end
% % % %     end
% % % %     D(ix) = D(ix)/nobjs(ix);        % to get the average
% % % % end
% % % % save
% % % % %% FUSE
% % % % %%% We have
% % % % %%% SC: view-rotation-aligned objects collected in bins, 
% % % % %%% SCpo: view-rotation-aligned and parameterization and object-space aligned collected in bins
% % % % %%% SA  : Angles that bring the objects from object-aligned back to original rotation around their own axis
% % % % %%% cc  : correspondence 
% % % % %%% D   : bin quality scores
% % % % load
% % % % T_av = [];
% % % % T = [];
% % % % for ix = 1:nbins,
% % % %     
% % % % %     tt = SC_ri{ix};
% % % % %     v1 = tt(1,:);v2 = tt(2,:);v3 = tt(3,:);
% % % %     
% % % %     vv = sum(SC_ri{ix})/size(SC_ri{ix},1);
% % % %     T_av = [T_av;vv];
% % % %     ang = SA{ix};   % take the first set of angles
% % % %     T(ix,:) = restore_configuration_from_r_inv(T_av(ix,:),ang(3,:));
% % % % end
% % % % save























