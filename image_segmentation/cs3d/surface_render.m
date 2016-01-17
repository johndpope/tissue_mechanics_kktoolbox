function [all_X, all_F, num_ren, Vol] = surface_render(segmented,V_max, V_min)
% %% Important: the input matrix Iw is assumed to be segmented (i.e. label matrix).
% %% The function returns a set of SHP coefficients for each segmented object.
% %% USAGE Example: [T, V_vec] = starting_shapes_local(Isegmented,I, 3, 100,500, 50);
% %% Notes: Needs  improvement in simplifying the surface triangulation prior
% %% to parameterization.
% %% Author: Khaled Khairy. Copyright EMBL-Heidelberg 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;

%% get the list of volumes and bounding boxes
stats = regionprops(segmented,'area');
num = length(stats);

for ix = 1:num,
    volume(ix) = stats(ix).Area;
end
volume = volume(:);
%%
one_step_render = 0;
num_ren = 0;    % number of rendered objects
all_X = {};
all_F = {};
Vol = [];
disp(['Number of labels: ' num2str(num)]);
thresh = 0.5;
if one_step_render
    %%% render in one step
    disp('Fast (but imprecise) rendering. Check ''surface_render.m'' for changing this.');
    num_ren = 1;
    segmented(segmented>0) = 1;
    data = smooth3(segmented, 'gaussian', [ 1 1 1]); % this is not always a good idea
    [all_F,all_X] = isosurface(data,thresh, 'noshare');       % generate the surface
    plot_mesh(all_X,all_F);
else
    %%% render shapes individually
    % [m,n,p] = size(segmented);  [x,y,z] = meshgrid(1:n,1:m,1:p);x = single(x); y = single(y); z = single(z); % use to initialize marching cubes
    
    BW = zeros(size(segmented));
    se = strel('ball',15,15);
    counter = 0;
    for t = 1:num
        if volume(t)<V_max && volume(t)>V_min,
            counter = counter + 1;
            disp(['Processing object: ' num2str(counter)]);
            BW = double(segmented==t);BW(BW>0) = 1;   % let's isolate this object (which could be present as some blob together with some smaller isolated pieces)
            %%% handle the case of noisy objects with many little pieces
            L = bwlabeln(BW, 6);                   % let us resegment this group
            while max(L(:))>1,     % i.e. if we have an object with many little pieces
                disp('Removing fragments');
                lp = 0;         % volume of largest piece
                ixp = 1;        % index of largest piece
                for ix = 1:max(L(:)),       % loop over these smaller pieces
                    if sum(L(:)==ix)>lp,
                        lp = sum(L(:)==ix);
                        ixp = ix;
                    end
                end
                BW(:) = 0;
                BW(L==ixp) = 1;
                %%% BW probably still has little handles. Let us close all gaps
%                BW = imclose(BW,se);
                BW = mat2gray(BW);
                BW(BW>0) = 1;
                Vol(counter) = sum(L(:)==ixp);
                L = bwlabeln(BW, 6);                   % let us resegment this group
            end
            disp(['Rendering object : ' num2str(t)]);
             [F,X] = isosurface(BW,thresh);       % generate the surface
%            [F,X] = kk_MarchingCubes(x,y,z,BW,thresh);
            all_X{counter} = X;
            all_F{counter} = F;
        end
    end
    %     if verbose,
    %         if ~isdeployed
    %             h = figure('visible','on');
    %             plot_all_X(all_X, all_F);
    %             saveas(h,'intermediate_result_rendering_object_surfaces.fig');
    %             print -dtiff -r600 intermediate_result_rendering_object_surfaces.tif;
    %             %close(h);
    %         end
    %     end
end

%% % delete objects with the wrong volume
indx = [];
for ix = 1:length(all_X)
    if ~(Vol(ix)<V_max && Vol(ix)>V_min),
        indx = [indx ix];
    end
end
all_X(indx) = [];
all_F(indx) = [];
Vol(indx) = [];
num_ren = length(all_F);









