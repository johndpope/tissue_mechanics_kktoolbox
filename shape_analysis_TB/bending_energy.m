function [Eb, A, V] = bending_energy(X_o)
% Calculate the bending energy of the object described by the clks
% initialize for speed by calling bending_energy_prep.m first

persistent Y_LKc HidAi dAi HidAi_mx dAij u v w crossqpr C
persistent V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix %LVE_ix_inv

%% precalculation of quantities needed for the fast Eb calculation%
L_max = get_L_max(X_o);
if isempty(Y_LKc) || size(Y_LKc,2)~=(L_max +1).^2,
    disp(' One-time bending energy calculation initialization');
    [X,C]=BuildSphere(3);
    [t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
    [x y z] = kk_sph2cart(t,p,1);
    Xc = double([x(:) y(:) z(:)]);
    Y_LKc = get_basis(t',p',length(t),get_L_max(X_o));
    
    [E,L] = edge_info(Xc,C);
    nvert = length(Xc);nfaces = length(C); nedges = nvert + nfaces - 2;
    LVE = sparse(nvert,nedges); VT = sparse(nvert,nfaces);V_e = sparse(nvert,nedges); VE_tr1 = sparse(nvert,nedges);
    VE_tr2 = sparse(nvert,nedges);V_far = sparse(nvert,nedges);dAij = sparse(nvert,nfaces);HidAi = sparse(nvert,1);
    dAi = sparse(nvert,1);HidAi_mx = sparse(nvert,nedges);u = sparse(nvert,1);v = sparse(nvert,1);w = sparse(nvert,1);
    crossqpr = sparse(nfaces,3);
    E = sort(E,2);      % sort the 2-vector according to its second dimension
    E = unique(E,'rows');% this should reduce the size of E by 50%.
    if nedges~=length(E),error('Number of edges is not correct');end
    %%% Now construct the matrix VE (i.e. vertex vs. edge connectivity)
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
        %%% which triangles does  the vertex ix belong to?
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
    
end
%% %
[xclks yclks zclks] = get_xyz_clks(X_o);
X = Y_LKc*[xclks(:) yclks(:) zclks(:)];
%%% calculation of the area and volume
u(:) = X(:,1); v(:) = X(:,2); w(:) = X(:,3);
crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2)); A = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;

% %%% Calculate H, dA and M (as in F.J. thesis): This requires the the sum of lengths of edges (j) for a particular vertex (i); Lij
HidAi_mx(LVE_ix) = sqrt((u(VE_ix)-u(V_e_ix)).^2 + (v(VE_ix)-v(V_e_ix)).^2+ (w(VE_ix)-w(V_e_ix)).^2).*...   %Lij
    real(acos(n(tr1_ix,1).*n(tr2_ix,1) + n(tr1_ix,2).*n(tr2_ix,2)+ n(tr1_ix,3).*n(tr2_ix,3))).*...    %thetaij
    (sign(n(tr1_ix,1).*u(V_far_ix) + n(tr1_ix,2).*v(V_far_ix) + n(tr1_ix,3).*w(V_far_ix) - ...
    n(tr1_ix,1).*u(VE_ix) - n(tr1_ix,2).*v(VE_ix) - n(tr1_ix,3).*w(VE_ix)));
HidAi(:) = 1/4 .* sum(HidAi_mx,2);dAij(LVT_ix) = F_areas(VT_ix)/3;dAi(:) = sum(dAij,2);M = sum(HidAi); H = HidAi./dAi;
dA  = 2 *sum(M);             % instantaneous area difference (assuming D = 1);


%%% Bending energy
Eb  = full(2* sum(HidAi.^2./dAi)); %bending energy as in F.J. thesis

Eb = Eb/8/pi;   % normalize to 1 ?