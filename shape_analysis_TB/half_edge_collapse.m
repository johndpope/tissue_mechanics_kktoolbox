function [X2, F2,L2, E, R, BC] = half_edge_collapse(X,F,L, niter)
%%% Do the half edge collapse as in Quicken et al.
%%% Input:  X = [x y z]      The actual coordinates of all points
%%%         F = The faces indices [vert1 vert 2 vert3] into X that give a
%%%         triangle. current
%%%         L = cell array of neighbors of each vertex. current
%%%
%%% Output:     X2: The new surface coordinate vector
%%%             F2: new faces after vertex removal
%%%             L2: new cell array of neighbors (corresponds to F2)
%%%             E : list of all edges (good for forces calculation)
%%%             R : list of removed indices
%%%             BC: barycentric coordinates of indices removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4,niter = 1;end; if size(F,1)==3, F = F';end;if size(X,1)==3, X = X';end;
X2 = X;F2 = F;L2 = L;verbose = 0;nvert_before =length(X);nfaces_before = length(F);
for iter = 1:niter,
    R = [];BC = [];nvert = length(X2);M = [];
    for vert = 1:nvert
        if verbose, if mod(vert,100)==0;disp(vert);end;end
        neighbors = L2{vert};    % get the list of neighbor indices
        fmemb = ismember(F2, vert);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
        if any(M == vert), if verbose, disp([num2str(vert) ' Vertex may not be removed: Neighbor of collapsed vertex']);end   % do not remove this vertex if it is on the list (M)
        elseif length(neighbors) <=3, if verbose,disp([num2str(vert) ' Vertex may not be removed: too few neighbors']); end;% only remove vertices that have at least 3 neighbors
        elseif isempty(fmemb), if verbose, disp([num2str(vert) ' Vertex has already been removed.']);end   % do not remove this vertex if it is on the list (M)
        else
            [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix

            cround = zeros(length(neighbors),1);    % this is where cround will be stored
            for ix = 1:length(neighbors)        % calculate cround for each possiblity of collapsing
                %     if length(L{neighbors(ix)}) >3,   % only consider neighbors with more than 3 neighbors
                C = F(ig,:);    %% the faces around vert
                target = neighbors(ix);     % this is the vertex that vert will collapse with.
                C(C==vert) = target;        % exchange all occurences of vert with target
                cc = sum(C==target,2);      % which triangles have two occurences of target?
                C = C(cc~=2,:);             % only keep the triangles that remain after half-edge collapse
                %         disp(C);
                for trindx = 1:size(C,1),   % loop over the triangles
                    ai = X2(C(trindx,2),:)-X2(C(trindx,1),:);bi = X2(C(trindx,3),:)-X2(C(trindx,1),:);ci = X2(C(trindx,3),:)-X2(C(trindx,2),:);
                    %         disp([ai;bi;ci]);
                    nai = norm(ai);nbi = norm(bi);nci = norm(ci);
                    s = 1/2* (nai + nbi + nci);
                    cround(ix) = cround(ix) + max([nai nbi nci]).*sqrt(s)./sqrt((s-nai)*(s-nbi)*(s-nci));
                end
                %     end
            end
            %%% now let us do the above for the collapse that minimizes cround
            target = neighbors(find(cround==min(cround)));     % this is the vertex that vert will collapse with.
            F2(F2==vert) = target(1);        % exchange all occurences of vert with target
            cc = sum(F2==target(1),2);      % which triangles have two occurences of target?
            F2 = F2(cc~=2,:);             % only keep the triangles that remain after half-edge collapse
            M = unique([M neighbors]);  % mark the neighbors as not removable

            if verbose,disp([num2str(vert) ' Removing vertex: collapsed with ' num2str(target(1))]);disp(cround');end;
            R = [R vert];
        end
    end
    [X2, F2]= remove_unused_vertices(X2,F2,R);
    [E,L2] = edge_info(X2,F2);
    nvert_after =length(X2);nfaces_after = length(F2);
    disp(['Before/after:  Vertices: ' num2str(nvert_before) '/' num2str(nvert_after)...
          '    Faces: ' num2str(nfaces_before) '/' num2str(nfaces_after)]);
end
%%%%%%%%%%%%%%%%%

