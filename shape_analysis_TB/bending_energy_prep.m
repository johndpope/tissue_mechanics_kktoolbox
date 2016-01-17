function bending_energy_prep(Xc, C)
%%% prepare quantities for bending energy calculation. To be used in
%%% conjunction with bending_energy.m.


%%%%%%%%%% Precalculate necessary quantities for the fast object curvature calculation
global HidAi dAi HidAi_mx dAij V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix u v w crossqpr

[E,L] = edge_info(Xc,C);
nvert = length(Xc);nfaces = length(C); nedges = nvert + nfaces - 2;
LVE = sparse(nvert,nedges); VT = sparse(nvert,nfaces);V_e = sparse(nvert,nedges); VE_tr1 = sparse(nvert,nedges);
VE_tr2 = sparse(nvert,nedges);V_far = sparse(nvert,nedges);dAij = sparse(nvert,nfaces);HidAi = sparse(nvert,1);
dAi = sparse(nvert,1);HidAi_mx = sparse(nvert,nedges);u = sparse(nvert,1);v = sparse(nvert,1);w = sparse(nvert,1);
crossqpr = sparse(nfaces,3);
E = sort(E,2);      % sort the 2-vector according to its second dimension
E = unique(E,'rows');% this should reduce the size of E by 50%.
if nedges~=length(E),error('Number of edges is not correct');end
%% Now construct the matrix VE (i.e. vertex vs. edge connectivity)
for ix = 1:nvert,   % loop over the vertices
    ememb = ismember(E, ix);% find the edges that vertex ix belongs to% find the rows in E where ix occurs
    [ig, jg] = ind2sub(size(ememb), find(ememb==1)); % ig is a column vector that indicates in which row of E we can find ix
    for ik = 1:length(ig),
        LVE(ix,ig(ik)) = ix;
        edg = E(ig(ik),:);      % 2-vector of the two vertices.
        V_e(ix,ig(ik)) = [edg(edg~=ix)];  % assigns the index of the other vertex to the matrix element
        g1 = ismember(C,edg(1));g2 = ismember(C,edg(2));tr = find(sum([g1 g2],2)==2);
        if length(tr)~=2,disp(tr);error('Could not determine unique two triangles of an edge.');end
        VE_tr1(ix,ig(ik)) = tr(1);
        VE_tr2(ix,ig(ik)) = tr(2); 
        tr2 = C(tr(2),:); edg = edg(edg~=ix);
        V_far(ix,ig(ik)) = [tr2(tr2~=ix & tr2~=edg)];
    end
    %% which triangles does  the vertex ix belong to?
    [r c] = ind2sub(size(C),find(C==ix)); % the rows define the triangles indexed into C (which defines the triangles)
    for ik = 1:length(r), VT(ix,r) = r;  end
end;
VE_ix       = (LVE(LVE~=0));
% warning off;LEV = (LVE); warning on
LVE_ix      = (find(LVE));% constructs the index vector that can be used to populate a matrix exactly where the edges exist
V_e_ix      = (V_e(V_e~=0));
V_far_ix    = (V_far(V_far~=0));
tr1_ix      = (VE_tr1(VE_tr1~=0)); % indices of triangle 1 positioned in VE
tr2_ix      = (VE_tr2(VE_tr1~=0));
VT_ix = VT(find(VT~=0));
LVT_ix = (find(VT~=0));
% disp(' Bending energy calculation initialized');
