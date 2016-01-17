function [all_X_out, all_F_out] = clean_all_XF(all_X, all_F, laplace_iter)
%%%%
%% clean and subdivide the shapes defined by the cell arrays all_X and
%% all_F. The shapes returned should be suitable for SHP fitting.

global ISO2MESH_TEMP
all_X_out = {};
all_F_out = {};
ISO2MESH_TEMP = './';
verbose = 0;
discard = [];
minVert = 4;
for ix = 1:size(all_X,2),   % loop over the shapes
    if verbose, disp(' ');
        disp(' |-------------------------------------|');
        disp(['    Cleaning surface number: ' num2str(ix)]);
        disp(' |-------------------------------------|');
    end
    X = all_X{ix};
    if size(X,1)>minVert,        % a minimum size should be in there
        F = all_F{ix};
        %%%%%%%%%%%% simplify using CGAL (computational geometry algorithms library), and then clean
        %        [X, F] = meshresample(X,F,0.5);
        [X, F] = meshcheckrepair(X,F,'duplicated');
        [X, F] = meshcheckrepair(X,F,'isolated');
        if~isdeployed, [X, F] = meshcheckrepair(X,F,'deep', num2str(ix));end
        facecell= finddisconnsurf(F);        % find disconnected surfaces
        if size(facecell,2)>1,
            if verbose, disp(' ****** multiple surface fragments found, picking largest!');end;
            F      = maxsurf(facecell);         % extract the largest surface
            %%%%%%%%%%%% simplify using CGAL (again) and clean mesh
            [X, F] = meshcheckrepair(X,F,'duplicated');
            [X, F] = meshcheckrepair(X,F,'isolated');
            if~isdeployed,[X, F] = meshcheckrepair(X,F,'deep', num2str(ix));end
        end
        if (length(X)*2-4~=length(F)), 
            warning('KK: Mesh does not appear to represent a closed shape!');
            %%% use convhull to close the mesh
            F = convhulln(X);
            [X, F] = meshcheckrepair(X,F,'deep', num2str(ix));
            [X, F] = subdivide(X,F,2);
            if verbose, if (length(X)*2-4~=length(F)), warning('********* KK: Mesh does not appear to represent a closed shape!*******');end;end
        end %#ok<WNTAG>
        
        %%%%%%%%%%%% generate the cell array of links in order to Laplacian smooth
        L = {};
        for jx = 1:length(X),   % loop over the vertices
            fmemb = ismember(F, jx);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
            [ig dummy] = ind2sub(size(fmemb), find(fmemb==1)); %#ok<NASGU> % ig is a column vector that indicates in which row of F we can find ix
            links = [];
            for ik = 1:length(ig),
                links = [links F(ig(ik),:)]; %#ok<AGROW>
            end
            L{jx} = unique(links(links~=jx));   %#ok<AGROW> % only record the links that are not jx
        end
        [X] = laplacian_smooth(X,F, L, laplace_iter);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_X_out{ix} = X;  % store the smoothed cleaned mesh
        all_F_out{ix} = F;
        status = fclose('all');
    end
end
