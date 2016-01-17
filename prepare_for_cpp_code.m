function prepare_for_cpp_code()
%% export theta, phi, F of triangular mesh on sphere with comma for c++ code embedding
[X, C] = mesh_gen(77,1);
[MESHINFO] = triangulated_props_mx_gen(X, C);
%[res, TP] = tri_tri_self_intersect(X,C);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
dlmwrite('phi.out',p', ',');
dlmwrite('theta.out',t', ',');
dlmwrite('F0.out',C(:,1)', ',');
dlmwrite('F1.out',C(:,2)', ',');
dlmwrite('F2.out',C(:,3)', ',');
%% export mesh information
% fid = fopen('uni900.tri', 'w');
fid = fopen('uni5929.tri', 'w');
for(ix = 1:numel(MESHINFO))
    n = MESHINFO{ix}.nm;
    fprintf(fid, '%d\n', n);
    for (nix = 1:n)
        fprintf(fid, '%d\t', MESHINFO{ix}.r(nix));
    end
    fprintf(fid, '\n');
    for (nix = 1:n)
        fprintf(fid, '%d\t%d\t%d\t%d\n', MESHINFO{ix}.indx(nix), MESHINFO{ix}.r1(nix), MESHINFO{ix}.r2(nix), MESHINFO{ix}.r3(nix));
    end
end
fclose(fid);

%% export shapes to ascii format
s = shp_surface('bowling_pin');fn = 'bowling_pin.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('sheet');fn = 'sheet.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('discocyte');fn = 'discocyte.shp3';%r_inv(s);
export_ascii(s, fn);
s = shp_surface('stomatocyte_01');fn = 'stomatocyte_01.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('stomatocyte_02');fn = 'stomatocyte_02.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('stomatocyte_03');fn = 'stomatocyte_03.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('plectrum');fn = 'plectrum.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('dros_embryo_01');fn = 'dros_embryo_01.shp3';export_ascii(r_inv(s), fn);
s = shp_surface('disc');fn = 'disc.shp3';export_ascii(r_inv(s), fn);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% generates the basis functions individually
L = 70;gdim = 40;
Y = [];%zeros((L+1)*(L+1),4);
tag = {};
ix = 1;
tag{ix} = 'L1Km1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L1K0';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L1K1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L2Km2';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L2Km1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L2K0';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L2K1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L2K2';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3Km3';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3Km2';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3Km1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3K0';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3K1';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3K2';Y(ix+1,ix) = 1.0;ix = ix + 1;
tag{ix} = 'L3K3';Y(ix+1,ix) = 1.0;ix = ix + 1;
s = shp_surface(L,gdim);
for(ix = 1:numel(tag))
    sh = sh_surface(L,gdim);sh.xc = sh_surface.tr_xc(Y(:,ix), L);
    sf = {tag{ix}, sh};
    s.sf{ix} = sf;
end
export_ascii(s,'basis_function.shp3');















function [MESHINFO] = triangulated_props_mx_gen(X, C)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Calculating properties from triangulation');
if size(X,1) ==3, X = X';end
A = 0; V = 0; H = 0; h = 0; v = 0; F_areas = 0;
%% quick calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2);
twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
A = sum(twoA)/2;                    % this is the total area
F_areas = twoA/2;                   % this is the vector of face areas
n = crossqpr./twoA(:,ones(3,1));
V = abs(sum(1/3*(dot(n,[x1 y1 z1], 2).*twoA./2)));
Vo = 4/3*pi*(A/4/pi)^(3/2);
v = V/Vo;
% disp([A V v]);

%%%% Now calculate the local mean curvature H at a vertex
%% Remember: Gauss curvature is the angle defect at a vector
%%          Mean curvature is edge length x dihedral angle at edge
H = zeros(size(u));         % the local mean curvature
N = zeros(length(u),3);     % store the average normal vector at a vertex
%normals_working = [];
K = zeros(size(H));         % the Gaussian curvature
MESHINFO = cell(length(X),1);
for ix = 1:length(X),   %loop over the vertices
    V_a = X(ix,:);
    %%% To calculate the curvature at some vertex V_a  we need to find the triangles that it is part of.
    [r c] = ind2sub(size(C),find(C==ix)); % the rows define the triangles indexed into C
    %% Calculate the dihedral angles of all unique permutaions of the triangle pairs (planes)
    [I2 I1] = ind2sub([length(r) length(r)],find( tril(ones(length(r)),-1)~=0)); % all possible combinations of triangles
    %% I2 and I1 are to be understood as pairs of indices into r. A value in r is an index into a row in C, which in turn contains indices
    %% into the points (rows) in X. so let us loop over the permutations and calculate the curvature H
    H(ix) = 0;
    theta = [];
    r1v = [];r2v = [];r3v = [];indxv = [];
    for I = 1:length(I2)% each permutation selects a triangle pair
        r1 = r(I2(I));r2 = r(I1(I));  % the two rows(triangles) in C that are considered
        tr1 = C(r1,:);tr2 = C(r2,:); % the two triangles (they contain indices into rows of X)
        % two of these indices are identical (this is the edge, the length of which we need to calculate we know that one of the vertices that are on the edge is V_a
        % now we need to find the other vertex V_e
        tr1r = tr1((tr1~=ix));tr2r = tr2((tr2~=ix));
        rvrs = (length(tr1)-1):-1:1;
        indx = max(tr1r.*(tr1r==tr2r) + tr1r(rvrs).*(tr1r(rvrs)==tr2r));    %  is the index of the other vertex V_e
%         r1vec = [r1vec r1];r2vec = [r2vec r2];
        if indx >0,  %i.e. if they share an edge
            r1v = [r1v;r1];
            r2v = [r2v;r2];
            r3v = [r3v;(tr2(tr2~=ix & tr2 ~=indx))];
            indxv = [indxv;indx];
            
            V_e = X(indx,:);
            V_far = X(tr2(tr2~=ix & tr2 ~=indx),:);
            % the distance between the two vertices of the common edge
            Lij   = sqrt((V_a(1)-V_e(1))^2 +(V_a(2)-V_e(2))^2+(V_a(3)-V_e(3))^2);
            %                 lijvec = [lijvec Lij];
            %% the surface normals are
            n1 = n(r1,:);n2 = n(r2,:);
            %normals_working = [normals_working;n1(:)';n2(:)'];
            %                 n1xvec = [n1xvec n1(1)];n2xvec = [n2xvec n2(1)];
            theta(I) = acos(dot(n1,n2));
            %%% we need to determine if they are convex or nonconvex
            % first let us complete the Hessian normal form of the planes of the two intersecting triangles
            P1 = -(n1(1)*V_a(1) + n1(2)*V_a(2) + n1(3)*V_a(3));
            P2 = -(n2(1)*V_a(1) + n2(2)*V_a(2) + n2(3)*V_a(3));
            %%% Now we need to calculate whether the far point of tr2 lies
            %%% in the half-space of the normal direction (local H is -ve)
            %%% or on the anti-normal direction (local curvature is +ve)
            s = sign(dot(n1,V_far)+P1);
            H(ix) = H(ix) + Lij * real(theta(I))/4 * (s); % as in F.J.thesis
        end
    end
    MESHINFO{ix}.nm = length(r);
    MESHINFO{ix}.r = r;     % indices of the faces that the vertex is member of
    MESHINFO{ix}.r1 = r1v;
    MESHINFO{ix}.r2 = r2v;
    MESHINFO{ix}.r3 = r3v;
    MESHINFO{ix}.indx = indxv;
    
    %normals_working = unique(normals_working,'rows');
    %normals_working = sum(normals_working)/length(normals_working);
    %N(ix,:) =normals_working;
    N(ix,:) = sum(n(r,:))/numel(r);
    %normals_working = [];
    %     (theta*180/pi)';
    % The average area of triangles around the vertex V_a is given by
    dA(ix) = sum(F_areas(r))/3;
    if dA(ix)==0,M(ix) = 0;warning('Probable error in area calculation II');
    else
        H(ix) = H(ix)/dA(ix);
        M(ix) = H(ix).*dA(ix);  % this is correct (as in F.J.thesis)
    end
end


% %     if nargout>9, K = real(K)./dA(:);k_g =sum(K.*dA.^2').*length(K)./A;end
% % %    if nargout>9, K = real(K)./dA(:);k_g = sum(K.*dA(:)).*length(K);end
% %
% %     h = sum((-M))./A;
% %     Eb  = 2* sum(M.^2./dA);   %bending energy as in F.J. thesis
% %     r = sqrt(A/4/pi);
% %     dA  = 2 *sum(-M);
% %     dAo = 8 * pi * r;
% %     da  = dA/dAo;
% %     Ebo = 8*pi;
% %     wb  = Eb/Ebo;
% %     D = 0.003;  % inner leaflet sepparation in microns
% %     str = sprintf('Area: %.5f \nVolume: %.5f\nReduced Volume: %.5f\nTotal mean curvature:%.5f\nE(bending): %.5f\nArea Difference: %.5f\n',...
% %         A, V,v,h,Eb/8/pi,dA*D);
% %     disp(str);
% %     % disp([A V v h Eb dA da]);
% %
% %     da = -da;

function [X, F,x, y, z,  A, V, v, F_areas, h, H, Eb, da] = mesh_gen(dim, flag)
A = 12.5664;

if(flag),
    %%%% old approximate method
    P = partsphere(dim^2);
    x = reshape(P(1,:),dim,dim);
    y = reshape(P(2,:),dim,dim);
    z = reshape(P(3,:),dim,dim);
else
    % % %%%% using subdivisions of icosahedron
    if dim>6, dim = 6;end
    [X,F]=BuildSphere(dim);
    [t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
    [x y z] = kk_sph2cart(t,p,sqrt(A/4/pi));
end

X = double([x(:) y(:) z(:)]);


[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
whos X x
C = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
F = C;

[A, V, v, F_areas, h, H, Eb, da] = triangulated_props(X, F, 0);

x = X(:,1);
y = X(:,2);
z = X(:,3);
% else
%     error('Too many triangles');
% end

