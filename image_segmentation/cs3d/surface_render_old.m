function [X, F, num_ren] = surface_render(segmented,V_max, V_min)
%% Important: the input matrix Iw is assumed to be segmented (i.e. label matrix).
%% The function returns a set of SHP coefficients for each segmented object.
%% USAGE Example: [T, V_vec] = starting_shapes_local(Isegmented,I, 3, 100,500, 50);
%% Notes: Needs massive improvement in simplifying the surface triangulation prior
%% to parameterization.
%% Author: Khaled Khairy. Copyright EMBL-Heidelberg 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
one_step_render = 0;
num_ren = 0;    % number of rendered objects
disp(['Final number of segmented objects found: ' num2str(max(segmented(:)))]);
axis on;axis equal;
view(3);
xlim([1 size(segmented,2)]);
ylim([1 size(segmented,1)]);
zlim([1 size(segmented,3)]);
light;
lighting gouraud;hold on;
thres = 0.5;
if one_step_render
    %%% render in one step
    disp('Fast (but imprecise) rendering. Check ''surface_render.m'' for changing this.');
    num_ren = 1;
    segmented(segmented>0) = 1;
    data = smooth3(segmented, 'gaussian', [ 5 5 5]); % this is not always a good idea
    [F,X] = isosurface(data,thres, 'noshare');       % generate the surface
    patch('Vertices', X, 'Faces', F,'FaceColor', 'r','EdgeColor','none','FaceAlpha',1);drawnow;

else
    %%% render shapes individually
    V = [];
    for t = max(segmented(:)):-1: 576%1
        BW = (segmented==t);BW(BW>0) = 1;
        L = bwlabeln(BW);
        V = [V sum(L(:)==1)]; %#ok<AGROW>
        disp(['Rendering object : ' num2str(t) ' with volume ' num2str(V(end))]);
        if V(end)<V_max && V(end)>V_min, %#ok<BDSCI>
            num_ren = num_ren + 1;
            tmp_im = (L==1);
            tmp_im = padarray(tmp_im,[1 1 1]);
            [F,X] = isosurface(squeeze(tmp_im),thres, 'noshare');       % generate the surface\
            %%%% check the surface for holes, fill them and fit to SHP
            if  length(X)*2-4~=length(F),disp([length(X)*2-4 length(F) (length(X)*2-4)-length(F)]);
                disp('Shape is not closed');
            end

            a = rand(1,3);      % random RGB color specification
            patch('Vertices', X, 'Faces', F,'FaceColor', round(a),'EdgeColor','none','FaceAlpha',1);
            drawnow;
        else
            str = sprintf('Object volume = %.2f , \n outside the volume bounds of max %.2f to min %.2f',...
                V(end), V_max, V_min);
            disp(str);
        end
    end
    hold off

end











