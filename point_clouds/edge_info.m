function [E,L,face_memb,Ld, Li] = edge_info(X,F)
%% generates E, which stores the indices (into X) of edge elements
%% generates L, which is a cell array of neighboring vertex indices
%% face_memb, cell array of indices representing faces that a vertex is member of.
%% Ld is the direct links: only if quadmesh
%% Li is the non-direct links: only if quadmesh
Ld = [];
Li = [];
% disp('Generating necessary edge and linkage information');
nvert = length(X);
nfaces = length(F);
if size(F,1)==3 || size(F,1) ==4,  F = F';end
if size(X,1)==3, X = X';end
if (nvert*2-4~=nfaces), disp('Mesh does not represent a closed shape (for triangulations)');end
% disp([(nvert*2-4) nfaces (nvert*2-4~=nfaces)]);
E = zeros(nvert,2);counter = 0;L = {};face_memb = {};
% determine the links L and the member faces

for ix = 1:length(X),   % loop over the vertices
%     if ~mod(ix,2000),disp(ix);end; 
    fmemb = ismember(F, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
    face_memb{ix} = find(sum(fmemb,2));
    [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
    links = [];
    for ik = 1:length(ig),
        links = [links F(ig(ik),:)];
    end
    L{ix} = unique(links(links~=ix));   % only record the links that are not ix
    % create the list of all edges
    llinks = L{ix};
    for jx = 1:length(llinks),
        counter = counter + 1;
        E(counter,:) = [ix llinks(jx)];
    end;
end;
%% 
if size(F,2)==4,        % then we have a quad mesh and also need direct and indirect links
    for ix = 1:length(X)
        dlinks = [];
        ilinks = [];
        u = face_memb{ix};
        for ixf = 1:numel(u)    % loop over the faces and determine the link type
            pos = find(F(u(ixf),:)==ix);    % find the position of the index in that face
            if pos ==1,
                dlinks = [dlinks F(u(ixf),2) F(u(ixf),4)]; 
                ilinks = [ilinks F(u(ixf),3)];
            elseif pos==2,
                dlinks = [dlinks F(u(ixf),1) F(u(ixf),3)]; 
                ilinks = [ilinks F(u(ixf),4)];
            elseif pos==3,
                dlinks = [dlinks F(u(ixf),2) F(u(ixf),4)]; 
                ilinks = [ilinks F(u(ixf),1)];
            elseif pos==4,
                dlinks = [dlinks F(u(ixf),3) F(u(ixf),1)]; 
                ilinks = [ilinks F(u(ixf),2)];
            end
        end
        dlinks = unique(dlinks);
        ilinks = unique(ilinks);
        Ld{ix} = dlinks;
        Li{ix} = ilinks;
    end
end

            
            
            
            
            
            
            
            
            
            
            
