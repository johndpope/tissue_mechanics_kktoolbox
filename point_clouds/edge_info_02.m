function [E,L,face_memb] = edge_info_02(X,F)
%% generates E, which stores the indices (into X) of edge elements
%% generates L, which is a cell array of neighboring vertex indices
%% face_memb, cell array of indices representing faces that a vertex is member of.
tr = TriRep(F,X(:,1), X(:,2), X(:,3));
% disp('Generating necessary edge and linkage information');
nvert = length(X);
nfaces = length(F);
if size(F,1)==3, F = F';end
if size(X,1)==3, X = X';end
if (nvert*2-4~=nfaces), disp('Mesh does not represent a closed shape (for triangulations)');end

face_memb = tr.vertexAttachments;
E = tr.edges;

for ix = 1:length(X),   % loop over the vertices
%     if ~mod(ix,2000),disp(ix);end; 
    fmemb =   face_memb{ix};
    links = [];
    for ik = 1:length(fmemb),
        links = [links F(fmemb(ik),:)];
    end
    L{ix} = unique(links(links~=ix));   % only record the links that are not ix
end;
