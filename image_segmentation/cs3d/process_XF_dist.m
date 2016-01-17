function [XALL, FALL, Xix, Fix] = ...
        process_XF_dist(XALL, FALL, Xix,Fix, laplace_iter,start, last)
%% clean and subdivide the shapes defined by the cell arrays all_X and
%% all_F. The shapes returned should be suitable for SHP fitting.

%%% decode the input by reconstructing all_X and all_F
all_X = {};
all_F = {};
Xcount = 1;
Fcount = 1;
for ix = 1:length(Xix)
    all_X{ix} = XALL(Xcount:Xix(ix),:);
    all_F{ix} = FALL(Fcount:Fix(ix),:);
    Xcount = Xix(ix)+1;
    Fcount = Fix(ix)+1;
end
%%%%
discard = [];
for ix = start:last,   % loop over the shapes
    disp(' ');
    disp(' |-------------------------------------|');
    disp(['    Cleaning surface number: ' num2str(ix)]);
    disp(' |-------------------------------------|');
    X = all_X{ix};
    if size(X,1)>40,        % a minimum size should be in there
        F = all_F{ix};
        %%%%%%%%%%%% simplify using CGAL (computational geometry
        %%%%%%%%%%%% algorithms library), and then clean
        [X, F] = meshresample(X,F,0.5);
        [X, F] = meshcheckrepair(X,F,'duplicated');
        [X, F] = meshcheckrepair(X,F,'isolated');
        if isempty(X)   % this shape is probably not worth looking at
            discard = [discard ix];
        else
            [X, F] = meshcheckrepair(X,F,'deep');
            facecell= finddisconnsurf(F);        % find disconnected surfaces
            if size(facecell,2)>1,
                disp(' ****** multiple surface fragments found, picking largest!');
                F      = maxsurf(facecell);         % extract the largest surface
                %%%%%%%%%%%% simplify using CGAL (again) and clean mesh
                [X, F] = meshcheckrepair(X,F,'duplicated');
                [X, F] = meshcheckrepair(X,F,'isolated');
                [X, F] = meshcheckrepair(X,F,'deep');
            end
            if (length(X)*2-4~=length(F)), warning('Mesh does not appear to represent a closed shape!');end %#ok<WNTAG>

            %%%%%%%%%%%% generate the cell array of links in order to Laplacian smooth
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
            all_X{ix} = X;  % store the smoothed cleaned mesh
            all_F{ix} = F;
            status = fclose('all');
        end
    end
end
%% discard the ones that are useless (should be maybe one or two)
% all_X(discard) = [];
% all_F(discard) = [];
%%% since dist matlab is not able to handle cell arrays properly we will
%%% fill them into X and F 3-vectors and handle the indexing
X = [];
F = [];
XALL = [];
FALL = [];
Xix = zeros(size(start:last));
Fix = zeros(size(start:last));
for ix = start:last,
    X = all_X{ix};
    F = all_F{ix};
    Xix(ix) = size(X,1);
    Fix(ix) = size(F,1);
    XALL = [XALL;X];
    FALL = [FALL;F];
       
end























