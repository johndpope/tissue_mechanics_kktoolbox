function [X1, C1, X2, C2] = shp_cut_x(X_o, s)
%%% cut the shp surface at the x-position s

%% step 0: generate a basis
[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = shp_surface.get_L_max(X_o);
nicos = 5;      % number of icosahedron subdivisions
[X, C]=surface_mesh.sphere_mesh_gen(nicos);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
%%% generate the basis
[L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
Y  = zeros(N, M, 'single');
for S = 1:length(L),
    Y(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
x = Y*xclks;
X = double(Y * [xclks yclks zclks]);

X1 = X;
X2 = X;
C1 = C;
C2 = C;

indx = find(x<=s);

%% find the indices of the triangles that intersect the x = s plane
large = 1e20;
Xs = [s large large; s -large large;s 0 -large];Cs = [1 2 3];
TP = [ones(size(C,1),1) [1:size(C,1)]'];
[trix] = tri_tri_intersect_local(Xs,Cs, X, C, TP);  % returns indices of intersecting triangles;
% plot_mesh(X,C(trix,:));

%% % put all vertices occuring in these triangles onto the x = plane
Xp = [];
for ix = 1:length(trix)
    f = C(trix(ix),:);
    cmy = [X(f(1),2)+X(f(2),2) + X(f(3),2)]/3;
    cmz = [X(f(1),3)+X(f(2),3) + X(f(3),3)]/3;
    Xp = [Xp; s cmy cmz];
%     Xp = [Xp; s X(f(1),2) X(f(1),3); s X(f(2),2) X(f(2),3); s X(f(3),2) X(f(3),3)];
end
% [B ind] = sort(Xp(:,2));
% Xp = Xp(ind,:);
degree = 3;
span = 5;
Xp(:,2) = smooth(Xp(:,2),span,'sgolay',degree);
% % Xp(:,3) = smooth(Xp(:,3),span,'sgolay',degree);

Xp(:,2) = smooth(Xp(:,2),span,'moving');
Xp(:,3) = smooth(Xp(:,3),span,'moving');

kk_plot3(Xp);view(-90,0)












function [indx] = tri_tri_intersect_local(X1,C1, X2, C2, TP)
%%% test whether any of the triangles defined in the meshes X1 C1 intersect 
%%% with any triangles in X2 C2. Call once to obtain TP and then include TP
%%% in the function call for speedup.
%%% 
%%% Output: res = 0 for no intersection, 1 for intersection of at least one triangle pair
%%%   Author: Khaled Khairy --- based on Moeller code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;

%%% generate the triangle lists T1 and T2 corresponding to those pairs
nTP = length(TP);
T1 = zeros(nTP,3, 'uint16');
T2 = zeros(nTP,3, 'uint16');

T1= [C1(TP(:,1),1) C1(TP(:,1),2) C1(TP(:,1),3)];
T2= [C2(TP(:,2),1) C2(TP(:,2),2) C2(TP(:,2),3)];

%% generate the vertices (starting at this point we have to have the current X calculated beforehand)
V21 = [X2(T2(:,1),1) X2(T2(:,1),2) X2(T2(:,1),3)]; % coordinates of the first vertex of the T2 triangles
V22 = [X2(T2(:,2),1) X2(T2(:,2),2) X2(T2(:,2),3)];
V23 = [X2(T2(:,3),1) X2(T2(:,3),2) X2(T2(:,3),3)];

V11 = [X1(T1(:,1),1) X1(T1(:,1),2) X1(T1(:,1),3)]; % coordinates of the first vertex of the T1 triangles
V12 = [X1(T1(:,2),1) X1(T1(:,2),2) X1(T1(:,2),3)];
V13 = [X1(T1(:,3),1) X1(T1(:,3),2) X1(T1(:,3),3)];

%% calculate surface normals and distances to construct the plane equations for T2
%N2 = cross( (V22-V21), (V23-V21) , 2);  % triangle normal
N2 = kk_cross((V22-V21),(V23-V21));
%d2 = -dot(N2,V21,2);                      % to complete the plane equation N2 . X + d2 = 0
d2 = -kk_dot(N2,V21);
%% Calculate distances of vertices of T1 triangles to T2 planes
dV11 = kk_dot(N2,V11) + d2;
dV12 = kk_dot(N2,V12) + d2;
dV13 = kk_dot(N2,V13) + d2;

%% for every pair test whether all vertices are on the plane of the corresponding triangle or not
%%% The point is that if the vertices are not on the plane and all
%%% distances have the same sign, then there is no intersection
Dzero1 = uint8(dV11==0) + uint8(dV12==0) + uint8(dV13==0);  % zero means no point is on the plane
Dsign1 = abs((sign(dV11)) + (sign(dV12)) + (sign(dV13)))~=3;% zero means all signs are the same (i.e.no intersection)
% one means possible intersection

%% calculate surface normals and distances to construct the plane equations for T1
N1 = kk_cross( (V12-V11), (V13-V11));  % triangle normal
d1 = -kk_dot(N1,V11);                      % to complete the plane equation N2 . X + d2 = 0
%% Calculate distances of vertices of T2 triangles to T1 planes
dV21 = kk_dot(N1,V21) + d1;
dV22 = kk_dot(N1,V22) + d1;
dV23 = kk_dot(N1,V23) + d1;
Dzero2 = uint8(dV21==0) + uint8(dV22==0) + uint8(dV23==0);  % zero means no point is on the plane
Dsign2 = abs((sign(dV21)) + (sign(dV22)) + (sign(dV23)))~=3;    % test whether they are all the same sign

% % %%% talk to me
% % if any(Dzero1), disp('Dzero1: found triangle(s) fully coplanar to corresponding triangle');end
% % if any(Dzero2), disp('Dzero2: found triangle(s) fully coplanar to corresponding triangle');end


indx = find(and((Dsign1==1),Dsign2==1)); % get the indices of the possibly intersecting pairs
% % % %%% resolve possible triangle intersection
% % % if ~isempty(indx)
% % %     %disp('possible triangle intersection detected');
% % %     N1I = N1(indx,:);N2I = N2(indx,:);
% % %     DI = kk_cross(N1I,N2I);  % this is the line  along the plane of both triangles, i.e. it intersects two edges in each triangle of a pair
% % %     
% % %     V11I = V11(indx,:);V12I = V12(indx,:);V13I = V13(indx,:);
% % %     V21I = V21(indx,:);V22I = V22(indx,:);V23I = V23(indx,:);
% % %     
% % %     dV11I = dV11(indx,:);dV12I = dV12(indx,:);dV13I = dV13(indx,:);
% % %     dV21I = dV21(indx,:);dV22I = dV22(indx,:);dV23I = dV23(indx,:);
% % %     
% % %     %%% now we need to find the pair of points that is on one side and separate
% % %     %%% it from the one point that is on the other side. We will rearrange the
% % %     %%% vertices so that V1 and V2 are on one side and V3 is on the other. We
% % %     %%% will then calculate the intersection
% % %     
% % %     %%% let us do this using for loops first
% % %     for ix = 1:size(DI,1),       % loop over the triangle pairs that will be checked for intersection
% % %         D = DI(ix,:);
% % %         V0 = V11I(ix,:);    % coordinates of the first vertex of the first triangle in the ixth pair to be compared
% % %         V1 = V12I(ix,:);    % coordinates of the second vertex of the first triangle in the ixth pair to be compared
% % %         V2 = V13I(ix,:);    % coordinates of the third vertex of the first triangle in the ixth pair to be compared
% % %         
% % %         U0 = V21I(ix,:);
% % %         U1 = V22I(ix,:);
% % %         U2 = V23I(ix,:);
% % %         
% % %         dv0 = dV11I(ix,:); dv1 = dV12I(ix,:);dv2 = dV13I(ix,:);
% % %         dv0dv1 = dv0*dv1;dv0dv2 = dv0*dv2;
% % %         
% % %         du0 = dV21I(ix,:); du1 = dV22I(ix,:);du2 = dV23I(ix,:);
% % %         du0du1 = du0*du1;du0du2 = du0*du2;
% % %         %/* compute and index to the largest component of D */
% % %         index = find(abs(D) == max(abs(D)));
% % %         index = index(1);
% % %         %/* this is the simplified projection onto L*/
% % %         vp0=V0(index);
% % %         vp1=V1(index);
% % %         vp2=V2(index);
% % %         
% % %         up0=U0(index);
% % %         up1=U1(index);
% % %         up2=U2(index);
% % %         
% % %         %%%%%%
% % %         [isect10 isect11] = COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2);
% % %         %/* compute interval for triangle 2 */
% % %         [isect20 isect21] = COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2);
% % % 
% % %         [isect10, isect11] = kk_sort(isect10,isect11);
% % %         [isect20, isect21] = kk_sort(isect20,isect21);
% % %          if(isect11<isect20 || isect21<isect10), 
% % %              %res = 0; 
% % %              %disp('No intersection based on intervals');
% % %          else
% % %              res = res + 1;
% % %              %disp('Intersection detected based on interval')
% % %              if flag, break;end
% % %          end
% % %     end
% % % else
% % %            %res = 0;
% % %             %disp('No intersection based on preliminary test');
% % % end
% % % function [a , b] = kk_sort(a,b)
% % % if a>b, c = a; a = b;b = c;end
% % %     
% % % function [isect0, isect1] = COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2)
% % % if(D0D1>0.0),
% % %     
% % %     %/* here we know that D0D2<=0.0 */                   \
% % %     %/* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
% % %     [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
% % %     
% % % elseif(D0D2>0.0),
% % %     
% % %     %/* here we know that d0d1<=0.0 */
% % %     [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
% % %     
% % % elseif(D1*D2>0.0 || D0~=0.0),
% % %     
% % %     %/* here we know that d0d1<=0.0 or that D0!=0.0 */
% % %     [isect0 isect1] =ISECT(VV0,VV1,VV2,D0,D1,D2);
% % %     
% % % elseif(D1~=0.0),
% % %     
% % %     [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
% % %     
% % % elseif(D2~=0.0),
% % %     
% % %     [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
% % %     
% % % else
% % %     
% % %     %/* triangles are coplanar */
% % %     disp('triangles are coplanar. Test for intersection not implemented yet');
% % %     isect0 = nan;
% % %     isect1 = nan;
% % %     %return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
% % %     
% % % end
% % %   
% % % 
% % % function [isect0, isect1] = ISECT(VV0,VV1,VV2,D0,D1,D2)
% % %               isect0=VV0+(VV1-VV0)*D0/(D0-D1); 
% % %               isect1=VV0+(VV2-VV0)*D0/(D0-D2);
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % % % %% % generate slice at s
% % % % % %%% we need to find theta and phi at the slice
% % % % % %%% first find theta
% % % % % F = TriScatteredInterp(p,double(x),t, 'linear');
% % % % % incr = 0.01;
% % % % % ang = min(p):incr:max(p);
% % % % % ts = F(ang,s*ones(size(ang)));
% % % % % ts(isnan(ts)) = [];
% % % % % 
% % % % % 
% % % % % %%% second find phi
% % % % % 
% % % % % F = TriScatteredInterp(t,double(x),p, 'linear');
% % % % % ps = F(ts,s*ones(size(ts)));
% % % % % 
% % % % % %%% evaluate basis at ps and ts
% % % % % [L, K] = shp_surface.indices_gen(1:(L_max + 1)^2); 
% % % % % M = length(L);
% % % % % N = length(ts);
% % % % % Ys  = zeros(N, M, 'single');
% % % % % for S = 1:length(L),
% % % % %     Ys(:,S) = sh_basis.ylk_bosh(L(S),K(S),ps,ts)'; % uses bosh version
% % % % % end;
% % % % % 
% % % % % X = double(Ys * [xclks yclks zclks]);
% % % % % kk_plot3(X);




% % t = t(indx);
% % p = p(indx);
% % Y  = zeros(length(t), M, 'single');
% % for S = 1:length(L),
% %     Y(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
% % end;
% % X = double(Y * [xclks yclks zclks]);
% % F = TriScatteredInterp(X(:,1), X(:,2), X(:,3));
% % 
% % ti = -100:1:100;
% % [qx,qy] = meshgrid(ti,ti);
% % qz = F(qx,qy);
% % mesh(qx,qy,qz);
% % hold on;
% % plot3(X(:,1), X(:,2), X(:,3),'.');




% % % % % % % % %%% Step 2: exclude all faces that have any of the excluded vertices in them
% % % % % % C1 = C;
% % % % % % indc = [];%zeros(length(C),1,'uint32');
% % % % % % Xi = (indx);
% % % % % % for ix = 1:size(C1,1)   % loop over the faces
% % % % % %    f = C1(ix,:);        % f has three indices (into X)
% % % % % %    if sum(ismember(f,Xi)),
% % % % % %        indc = [indc ix];
% % % % % %    end
% % % % % % %    
% % % % % % %    t1 = sum(f(1)*ones(size(Xi))==Xi);
% % % % % % %    t2 = sum(f(2)*ones(size(Xi))==Xi);
% % % % % % %    t3 = sum(f(3)*ones(size(Xi))==Xi);
% % % % % % %    if sum([t1(:)' t2(:)' t3(:)'])>0,
% % % % % % %        disp('found');disp(ix);
% % % % % % %        indc(ix) = ix;
% % % % % % %    end  % mark this face for deletion
% % % % % % end
% % % % % % C1(indc,:) = [];
% % % % % % 
% % % % % % plot_mesh(X,C1);





