function [X_o Y_LK] = ade_ms_precalc(A_o, L_max, gdim, A_ms, v_red_ms, lambda_shear)
%
%

%%% Get starting values
disp(' Calculating initial values');
global C Y_LK HidAi dAi HidAi_mx dAij
global V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix u v w crossqpr VertTri
global Y_LK_C1 Y_LK_C2 Y_LK_C3 
%%%% prepare initial mesh
nico = 3;
[X, C,x, y, z,  A, V, v, F_areas, h, H, Eb, da] = mesh_gen(nico, A_o);nvert = length(X);
[t p Radius] = kk_cart2sph(x, y, z);
[A, V, h, v, xclks, yclks, zclks] = sh_projection(gdim, L_max, L_max, L_max, x', y', z', t', p', 0, 0);
X_o = [xclks(:)' yclks(:)' zclks(:)'];
Y_LK = get_basis(t',p',gdim,L_max); % just to get Y_LK;
    
counter = 0;for ix = 1:nvert,   % loop over the vertices
        fmemb = ismember(C, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
        [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
        links = [];
        for ik = 1:length(ig), 
            links = [links C(ig(ik),:)];
        end
        L{ix} = unique(links(links~=ix));   % only record the links that are not ix
        % create the list of all edges
        llinks = L{ix};
        for jx = 1:length(llinks),
            counter = counter + 1;
            E(counter,:) = [ix llinks(jx)];
        end;
end;
[L, K ] = indices_gen(1:(L_max + 1)^2); M = length(L);N = length(x);
Y_LK_tri  = zeros(N, M);for S = 1:length(L), Y_LK_tri(:,S) = ylk_cos_sin_old(L(S),K(S),p',t')';end;
%%%% For the vectorized area and volume calculation
Y_LK_C1  = zeros(length(C), M);for S = 1:length(L), Y_LK_C1(:,S) = ylk_cos_sin_old(L(S),K(S),p(C(:,1))',t(C(:,1))')';end;
Y_LK_C2  = zeros(length(C), M);for S = 1:length(L), Y_LK_C2(:,S) = ylk_cos_sin_old(L(S),K(S),p(C(:,2))',t(C(:,2))')';end;
Y_LK_C3  = zeros(length(C), M);for S = 1:length(L), Y_LK_C3(:,S) = ylk_cos_sin_old(L(S),K(S),p(C(:,3))',t(C(:,3))')';end;

%%%% precalculate necessary quantities
nvert
nfaces = length(C); nedges = nvert + nfaces - 2;
LVE = sparse(nvert,nedges); VT = sparse(nvert,nfaces);
V_e = sparse(nvert,nedges); VE_tr1 = sparse(nvert,nedges);
VE_tr2 = sparse(nvert,nedges);V_far = sparse(nvert,nedges);
dAij = sparse(nvert,nfaces);HidAi = sparse(nvert,1);
dAi = sparse(nvert,1);HidAi_mx = sparse(nvert,nedges);
u = sparse(nvert,1);v = sparse(nvert,1);w = sparse(nvert,1);
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
        if length(tr)~=2,error('Could not derermine unique two triangles of an edge.');end
        VE_tr1(ix,ig(ik)) = tr(1);
        VE_tr2(ix,ig(ik)) = tr(2); 
        tr2 = C(tr(2),:); edg = edg(edg~=ix);
        V_far(ix,ig(ik)) = [tr2(tr2~=ix & tr2~=edg)];
    end
    %% which triangles does  the vertex ix belong to?
    [r c] = ind2sub(size(C),find(C==ix)); % the rows define the triangles indexed into C (which defines the triangles)
    VertTri{ix} = r;
    for ik = 1:length(r)
        VT(ix,r) = r;
    end
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
X_o = [xclks; yclks; zclks]';
disp(' end of initialization of the ADE energy calculation');
clear nedges phi theta tr tr2 wp wt 
clear x y z xclks yclks zclks sing r llinks links jx jg ix ik
clear ig h gdim g2 g1 fmemb ememb edg counter c V_far 
clear V_e VT VE_tr1 VE_tr2 V LEV E A


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Precalculate quantities for the MS energy
if lambda_shear
global F_areas_o invDM DM 
%%% generate the ellipse
% v_target = 0.95; A_So = A_o;
Vo = 4/3*pi*(A_ms/4/pi).^(3/2);
V_target = v_red_ms * Vo;    % this is the actual volume of S_o in microns^3


%%%%%%%%%%%%%%% Oblate
[a, c, A_oblate, Vol, v_ell] = ellipsoid_dimensions_gen(A_ms, v_red_ms, 1);
if ~(a>=c), error('a must be larger than c');end
str  = sprintf('a:%.4f\nc:%.4f\nTarget : Actual Area: %.4f : %.4f\n Target : Actual volume:%.4f : %.4f\n Reduced volume: %.4f',...
    a,c,A_ms,A_oblate,V_target,Vol,v_ell);
disp(str);

% e = sqrt(1-(c^2/a^2));fac = log((1+e)/(1-e));% A_oblate = 2*pi*a^2 + pi * c^2/e*fac;  % mathworks% Vol = 4/3*pi*a^2 * c;
[x y z] = kk_sph2cart(t,p,1);
u = a * x; v = a * y; w = c * z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],...
    [u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));  F_areas_o = twoA/2;n = crossqpr./twoA(:,ones(3,1));
%%%%% build the necessary transformations as in Computer Graphics p. 217.
%%% Get the vertices 
V1 = [u(C(:,1)) v(C(:,1)) w(C(:,1))];V2 = [u(C(:,2)) v(C(:,2)) w(C(:,2))];V3 = [u(C(:,3)) v(C(:,3)) w(C(:,3))];
%%% construct the translation matrix that brings V1 to the origin
T  = speye(4*nfaces);
diagx1 = zeros(1,size(T,1));diagx1(4:4:end) = -V1(:,1);T = T + spdiags(diagx1',3,size(T,1), size(T,2));
diagy1 = zeros(1,size(T,1));diagy1(4:4:end) = -V1(:,2);T = T + spdiags(diagy1',2,size(T,1), size(T,2));
diagz1 = zeros(1,size(T,1));diagz1(4:4:end) = -V1(:,3);T = T + spdiags(diagz1',1,size(T,1), size(T,2));
%%% construct the necessary rotation matrix
R = speye(4*nfaces); p1p2 = V2-V1;p1p3 = V3-V1;
normp1p2 = sqrt(p1p2(:,1).^2 + p1p2(:,2).^2 + p1p2(:,3).^2);
rz = p1p2./normp1p2(:,ones(3,1));
crossp1p3p1p2 = cross(p1p3,p1p2, 2);normcrossp1p3p1p2 = sqrt(crossp1p3p1p2(:,1).^2 + crossp1p3p1p2(:,2).^2 + crossp1p3p1p2(:,3).^2);
rx = crossp1p3p1p2./normcrossp1p3p1p2(:,ones(3,1));ry = cross(rz, rx, 2);
diag_vec = ones(1,4 * nfaces);diag_vec(1:4:end) =rx(:,1);diag_vec(2:4:end) = ry(:,2);diag_vec(3:4:end) = rz(:,3);
R = R + spdiags(diag_vec',0,4*nfaces, 4*nfaces);
diag_vec_1 = zeros(1,4 * nfaces);diag_vec_1(2:4:end) = rx(:,2);diag_vec_1(3:4:end) = ry(:,3);R = R + spdiags(diag_vec_1',1,size(R,1), size(R,2));
diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_m1(1:4:end) = ry(:,1);diag_vec_m1(2:4:end) = rz(:,2);R = R + spdiags(diag_vec_m1',-1,size(R,1), size(R,2));
diag_vec_2 = zeros(1,4 * nfaces);diag_vec_2(3:4:end) = rx(:,3);R = R + spdiags(diag_vec_2',2,size(R,1), size(R,2));
diag_vec_m2 = zeros(1,4 * nfaces);diag_vec_m2(1:4:end) = rz(:,1);R = R + spdiags(diag_vec_m2',-2,size(R,1), size(R,2));
R = R - speye(4 * nfaces);
V1pr = reshape([V1 ones(length(V1),1)]',length(V1)*4,1);
V2pr = reshape([V2 ones(length(V2),1)]',length(V2)*4,1);
V3pr = reshape([V3 ones(length(V3),1)]',length(V3)*4,1);
V1r = R*T*V1pr; V2r = R * T *V2pr; V3r = R * T * V3pr; % transform
V1r = [reshape(V1r,4,length(V1))]';V2r = [reshape(V2r,4,length(V2))]';V3r = [reshape(V3r,4,length(V3))]';
%%% All vertices are rotated now, and all x-coordinates are zero.now we omit the x-coordinate and look at the triangles in 2-D only. note that all V1r are translated to zero.
V2r = [V2r(:,2) V2r(:,3)]; V3r = [V3r(:,2) V3r(:,3)];
%%% use the rotated vertices to build the necessary matrices for the calculation of the deformation
DM = speye(2*nfaces);%invDM = sparse(nfaces*2);
diag_vec = zeros(1,2*nfaces);diag_vec(1:2:end) =V2r(:,1);diag_vec(2:2:end) = V3r(:,2);DM = DM + spdiags(diag_vec',0,size(DM,1), size(DM,2));
diag_vec_1 = zeros(1,2*nfaces);diag_vec_1(2:2:end)=V3r(:,1);DM = DM + spdiags(diag_vec_1',1,size(DM,1), size(DM,2));
diag_vec_m1 = zeros(1,2*nfaces);diag_vec_m1(1:2:end)=V2r(:,2);DM = DM + spdiags(diag_vec_m1',-1,size(DM,1), size(DM,2));
DM = DM-speye(size(DM));

% invDM = inv(DM); 
%%% manually invert DM to make the inversion possible for fine meshes
indDM = 1:2:2*nfaces;
invDM = sparse(size(DM));
for ix = 1:nfaces,      % loop over the triangles (direction cosines rotation, and prep. DM)
    Dm = DM(indDM(ix):indDM(ix)+1,indDM(ix):indDM(ix)+1);
    invDM(indDM(ix):indDM(ix)+1,indDM(ix):indDM(ix)+1) = inv(Dm);
end


clear T R normp1p2 rz crossp1p3p1p2
clear  V1pr V2pr V3pr V1r V2r V3r
%%
global diagx1 diagy1 diagz1 SPI diag_vec_1 diag_vec diag_vec_m1 diag_vec_2 diag_vec_m2
global diag_vec_ds diag_vec_1_ds diag_vec_m1_ds SPI2
diagx1 = zeros(1,4*nfaces);diagy1 = zeros(1,4*nfaces);diagz1 = zeros(1,4*nfaces);
SPI = speye(4*nfaces); SPI2 = speye(2*nfaces);
diag_vec_1 = zeros(1,4 * nfaces);diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_2 = zeros(1,4 * nfaces);diag_vec_m2 = zeros(1,4 * nfaces);
diag_vec_ds = zeros(1,2*nfaces);diag_vec_1_ds = zeros(1,2*nfaces);diag_vec_m1_ds = zeros(1,2*nfaces);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = sparse(nvert,1);v = sparse(nvert,1);w = sparse(nvert,1);    % since these were used above in another way




save data_precalc