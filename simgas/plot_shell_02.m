function plot_shell_02(X,F, d, yval)
% the usual X and C, and the thickness d. val is the y value of the xz cutting plane

% % %% translate to center of mass (optional)
cm = sum(X,1)./size(X,1);
cm = cm(ones(1,size(X,1)),:);
X(:,1) = X(:,1)-cm(:,1);
X(:,3) = X(:,3)-cm(:,3);

% patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 0.8);  % the middle surface
%%% now that we have the normals at each point we calculate the surface that
%%% is a distance d removed from it
[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, F, 0);
X_in = X + N .* d/2;  %%% now generate the inner surface
X_out = X - N .* d/2;  %%% now generate the outer surface
%% find the triangles that intersect that surface
large = 1e3;
X2 = [-large yval -large; 0 yval large; large yval -large];F2 = [1 2 3];
[int_indx, TP] = tri_plane_intersect(X_in,F, X2, F2, [], 0); % determine the triangle indices that intersect the plane
%% of the intersecting triangles, project the vertex that is closest to the plane yval onto yval
contour = [];
for ix = 1:numel(int_indx),         % loop over intersecting triangles
   %%% find the vertex closest to the plane yval
   tri = F(int_indx(ix),:);
   dmin = large;
   minix = 1;
   for dix = 1:3
       v = X(tri(dix),:);d = abs(v(2)-yval);
       if d<dmin, dmin = d;minix = dix;end
   end
   %%% set the yvalue of that vertex to yval
   v = X(tri(minix),:);
   X(tri(minix),:) = [v(1) yval v(3)];
   %%% store this contour point as a 2-vector
   contour = [contour;v(1) v(3)];   % we only need to store the x and z values
end
% %plot(contour(:,1)', contour(:,2)','*');
%%  delete all triangles with vertices of values larger than yval, except those that only have one value
del_tri = [];
ylix_vec = [];
for ix = 1:size(F,1),   % loop over the triangles
    tri = F(ix,:);
    cum = 0;
    ylix = 0;
    for dix = 1:3
        v = X(tri(dix),:);
        if v(2)>yval, 
            cum = cum +1;
            ylix = dix;
        end
    end
    if cum>1, 
        del_tri = [del_tri;ix];
    elseif cum==1
       %%% the vertex with yval larger is stored in ylix
       ylix_vec = [ylix_vec;tri(ylix)]; % store the index of the point whose y value will be projected
    end
end
F(del_tri,:) = [];
%% now project the one vertex in the triangles identified in the last step 
for ix = 1:numel(ylix_vec)
       v1 = X(ylix_vec(ix),:);
       X(ylix_vec(ix),:) = [v1(1) yval v1(3)];
end
%%
plot_mesh(X,F);
% %% plot
% patch('vertices',X_in,'faces',F,'FaceColor','blue', 'FaceAlpha', 0.8);
% patch('vertices',X_out,'faces',F,'FaceColor','blue', 'FaceAlpha', 0.8);

%% find the points that are 

%%
function [int_indx, TP] = tri_plane_intersect(X1,C1, X2, C2, TP, flag)
%%% find which of the triangles defined in the meshes X1 C1 intersect
%%% with any triangles in X2 C2. Call once to obtain TP and then include TP
%%% in the function call for speedup.
%%%
%%% Output: res = 0 for no intersection, 1 for intersection of at least one triangle pair
%%%   Author: Khaled Khairy --- based on Moeller code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;
int_indx = [];
if nargin<6, flag = 0;end
if nargin < 5 || isempty(TP)
    %% generate triangle pair list
    %%% although the generation of this is slow it needs to be done only once
    %%% for other calls include TP in the call.
    TP = [];
    disp('Building triangle pair database....')
    for ix = 1:size(C1,1)
        %disp([num2str(ix) ' of ' num2str(size(C1,1))]);
        TP = [TP;[ix 1]];
    end
    disp('Done!');
end
%% % generate the triangle lists T1 and T2 corresponding to those pairs
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
%% % resolve possible triangle intersection
if ~isempty(indx)
    %disp('possible triangle intersection detected');
    N1I = N1(indx,:);N2I = N2(indx,:);
    DI = kk_cross(N1I,N2I);  % this is the line  along the plane of both triangles, i.e. it intersects two edges in each triangle of a pair
    
    V11I = V11(indx,:);V12I = V12(indx,:);V13I = V13(indx,:);
    V21I = V21(indx,:);V22I = V22(indx,:);V23I = V23(indx,:);
    
    dV11I = dV11(indx,:);dV12I = dV12(indx,:);dV13I = dV13(indx,:);
    dV21I = dV21(indx,:);dV22I = dV22(indx,:);dV23I = dV23(indx,:);
    
    %%% now we need to find the pair of points that is on one side and separate
    %%% it from the one point that is on the other side. We will rearrange the
    %%% vertices so that V1 and V2 are on one side and V3 is on the other. We
    %%% will then calculate the intersection
    
    %%% let us do this using for loops first
    for ix = 1:size(DI,1),       % loop over the triangle pairs that will be checked for intersection
        D = DI(ix,:);
        V0 = V11I(ix,:);    % coordinates of the first vertex of the first triangle in the ixth pair to be compared
        V1 = V12I(ix,:);    % coordinates of the second vertex of the first triangle in the ixth pair to be compared
        V2 = V13I(ix,:);    % coordinates of the third vertex of the first triangle in the ixth pair to be compared
        
        U0 = V21I(ix,:);
        U1 = V22I(ix,:);
        U2 = V23I(ix,:);
        
        dv0 = dV11I(ix,:); dv1 = dV12I(ix,:);dv2 = dV13I(ix,:);
        dv0dv1 = dv0*dv1;dv0dv2 = dv0*dv2;
        
        du0 = dV21I(ix,:); du1 = dV22I(ix,:);du2 = dV23I(ix,:);
        du0du1 = du0*du1;du0du2 = du0*du2;
        %/* compute and index to the largest component of D */
        index = find(abs(D) == max(abs(D)));
        index = index(1);
        %/* this is the simplified projection onto L*/
        vp0=V0(index);
        vp1=V1(index);
        vp2=V2(index);
        
        up0=U0(index);
        up1=U1(index);
        up2=U2(index);
        
        %%%%%%
        [isect10 isect11] = COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2);
        %/* compute interval for triangle 2 */
        [isect20 isect21] = COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2);
        
        [isect10, isect11] = kk_sort(isect10,isect11);
        [isect20, isect21] = kk_sort(isect20,isect21);
        if(isect11<isect20 || isect21<isect10),
            %res = 0;
            %disp('No intersection based on intervals');
        else
            res = res + 1;
            int_indx = [int_indx;indx(ix)];
            %disp('Intersection detected based on interval')
            if flag, break;end
        end
    end
else
    %res = 0;
    %disp('No intersection based on preliminary test');
end
function [a , b] = kk_sort(a,b)
if a>b, c = a; a = b;b = c;end
function [isect0, isect1] = COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2)
if(D0D1>0.0),
    
    %/* here we know that D0D2<=0.0 */                   \
    %/* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
    
elseif(D0D2>0.0),
    
    %/* here we know that d0d1<=0.0 */
    [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
    
elseif(D1*D2>0.0 || D0~=0.0),
    
    %/* here we know that d0d1<=0.0 or that D0!=0.0 */
    [isect0 isect1] =ISECT(VV0,VV1,VV2,D0,D1,D2);
    
elseif(D1~=0.0),
    
    [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
    
elseif(D2~=0.0),
    
    [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
    
else
    
    %/* triangles are coplanar */
    disp('triangles are coplanar. Test for intersection not implemented yet');
    isect0 = nan;
    isect1 = nan;
    %return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
    
end
function [isect0, isect1] = ISECT(VV0,VV1,VV2,D0,D1,D2)
isect0=VV0+(VV1-VV0)*D0/(D0-D1);
isect1=VV0+(VV2-VV0)*D0/(D0-D2);
%% The C-code from Moeller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% MOELLER code
% % /* Triangle/triangle intersection test routine,
% %  * by Tomas Moller, 1997.
% %  * See article "A Fast Triangle-Triangle Intersection Test",
% %  * Journal of Graphics Tools, 2(2), 1997
% %  * updated: 2001-06-20 (added line of intersection)
% %  *
% %  * int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
% %  *                       float U0[3],float U1[3],float U2[3])
% %  *
% %  * parameters: vertices of triangle 1: V0,V1,V2
% %  *             vertices of triangle 2: U0,U1,U2
% %  * result    : returns 1 if the triangles intersect, otherwise 0
% %  *
% %  * Here is a version withouts divisions (a little faster)
% %  * int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
% %  *                      float U0[3],float U1[3],float U2[3]);
% %  *
% %  * This version computes the line of intersection as well (if they are not coplanar):
% %  * int tri_tri_intersect_with_isectline(float V0[3],float V1[3],float V2[3],
% %  *				        float U0[3],float U1[3],float U2[3],int *coplanar,
% %  *				        float isectpt1[3],float isectpt2[3]);
% %  * coplanar returns whether the tris are coplanar
% %  * isectpt1, isectpt2 are the endpoints of the line of intersection
% %  */
% %
% % #include <math.h>
% %
% % #define FABS(x) ((float)fabs(x))        /* implement as is fastest on your machine */
% %
% % /* if USE_EPSILON_TEST is true then we do a check:
% %          if |dv|<EPSILON then dv=0.0;
% %    else no check is done (which is less robust)
% % */
% % #define USE_EPSILON_TEST TRUE
% % #define EPSILON 0.000001
% %
% %
% % /* some macros */
% % #define CROSS(dest,v1,v2)                      \
% %               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
% %               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
% %               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
% %
% % #define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
% %
% % #define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2];
% %
% % #define ADD(dest,v1,v2) dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2];
% %
% % #define MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2];
% %
% % #define SET(dest,src) dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2];
% %
% % /* sort so that a<=b */
% % #define SORT(a,b)       \
% %              if(a>b)    \
% %              {          \
% %                float c; \
% %                c=a;     \
% %                a=b;     \
% %                b=c;     \
% %              }
% %
% % #define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
% %               isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
% %               isect1=VV0+(VV2-VV0)*D0/(D0-D2);
% %
% %
% % #define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
% %   if(D0D1>0.0f)                                         \
% %   {                                                     \
% %     /* here we know that D0D2<=0.0 */                   \
% %     /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
% %     ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
% %   }                                                     \
% %   else if(D0D2>0.0f)                                    \
% %   {                                                     \
% %     /* here we know that d0d1<=0.0 */                   \
% %     ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
% %   }                                                     \
% %   else if(D1*D2>0.0f || D0!=0.0f)                       \
% %   {                                                     \
% %     /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
% %     ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
% %   }                                                     \
% %   else if(D1!=0.0f)                                     \
% %   {                                                     \
% %     ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
% %   }                                                     \
% %   else if(D2!=0.0f)                                     \
% %   {                                                     \
% %     ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
% %   }                                                     \
% %   else                                                  \
% %   {                                                     \
% %     /* triangles are coplanar */                        \
% %     return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
% %   }
% %
% %
% %
% % /* this edge to edge test is based on Franlin Antonio's gem:
% %    "Faster Line Segment Intersection", in Graphics Gems III,
% %    pp. 199-202 */
% % #define EDGE_EDGE_TEST(V0,U0,U1)                      \
% %   Bx=U0[i0]-U1[i0];                                   \
% %   By=U0[i1]-U1[i1];                                   \
% %   Cx=V0[i0]-U0[i0];                                   \
% %   Cy=V0[i1]-U0[i1];                                   \
% %   f=Ay*Bx-Ax*By;                                      \
% %   d=By*Cx-Bx*Cy;                                      \
% %   if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
% %   {                                                   \
% %     e=Ax*Cy-Ay*Cx;                                    \
% %     if(f>0)                                           \
% %     {                                                 \
% %       if(e>=0 && e<=f) return 1;                      \
% %     }                                                 \
% %     else                                              \
% %     {                                                 \
% %       if(e<=0 && e>=f) return 1;                      \
% %     }                                                 \
% %   }
% %
% % #define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
% % {                                              \
% %   float Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
% %   Ax=V1[i0]-V0[i0];                            \
% %   Ay=V1[i1]-V0[i1];                            \
% %   /* test edge U0,U1 against V0,V1 */          \
% %   EDGE_EDGE_TEST(V0,U0,U1);                    \
% %   /* test edge U1,U2 against V0,V1 */          \
% %   EDGE_EDGE_TEST(V0,U1,U2);                    \
% %   /* test edge U2,U1 against V0,V1 */          \
% %   EDGE_EDGE_TEST(V0,U2,U0);                    \
% % }
% %
% % #define POINT_IN_TRI(V0,U0,U1,U2)           \
% % {                                           \
% %   float a,b,c,d0,d1,d2;                     \
% %   /* is T1 completly inside T2? */          \
% %   /* check if V0 is inside tri(U0,U1,U2) */ \
% %   a=U1[i1]-U0[i1];                          \
% %   b=-(U1[i0]-U0[i0]);                       \
% %   c=-a*U0[i0]-b*U0[i1];                     \
% %   d0=a*V0[i0]+b*V0[i1]+c;                   \
% %                                             \
% %   a=U2[i1]-U1[i1];                          \
% %   b=-(U2[i0]-U1[i0]);                       \
% %   c=-a*U1[i0]-b*U1[i1];                     \
% %   d1=a*V0[i0]+b*V0[i1]+c;                   \
% %                                             \
% %   a=U0[i1]-U2[i1];                          \
% %   b=-(U0[i0]-U2[i0]);                       \
% %   c=-a*U2[i0]-b*U2[i1];                     \
% %   d2=a*V0[i0]+b*V0[i1]+c;                   \
% %   if(d0*d1>0.0)                             \
% %   {                                         \
% %     if(d0*d2>0.0) return 1;                 \
% %   }                                         \
% % }
% %
% % int coplanar_tri_tri(float N[3],float V0[3],float V1[3],float V2[3],
% %                      float U0[3],float U1[3],float U2[3])
% % {
% %    float A[3];
% %    short i0,i1;
% %    /* first project onto an axis-aligned plane, that maximizes the area */
% %    /* of the triangles, compute indices: i0,i1. */
% %    A[0]=fabs(N[0]);
% %    A[1]=fabs(N[1]);
% %    A[2]=fabs(N[2]);
% %    if(A[0]>A[1])
% %    {
% %       if(A[0]>A[2])
% %       {
% %           i0=1;      /* A[0] is greatest */
% %           i1=2;
% %       }
% %       else
% %       {
% %           i0=0;      /* A[2] is greatest */
% %           i1=1;
% %       }
% %    }
% %    else   /* A[0]<=A[1] */
% %    {
% %       if(A[2]>A[1])
% %       {
% %           i0=0;      /* A[2] is greatest */
% %           i1=1;
% %       }
% %       else
% %       {
% %           i0=0;      /* A[1] is greatest */
% %           i1=2;
% %       }
% %     }
% %
% %     /* test all edges of triangle 1 against the edges of triangle 2 */
% %     EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
% %     EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
% %     EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
% %
% %     /* finally, test if tri1 is totally contained in tri2 or vice versa */
% %     POINT_IN_TRI(V0,U0,U1,U2);
% %     POINT_IN_TRI(U0,V0,V1,V2);
% %
% %     return 0;
% % }
% %
% %
% % int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
% %                       float U0[3],float U1[3],float U2[3])
% % {
% %   float E1[3],E2[3];
% %   float N1[3],N2[3],d1,d2;
% %   float du0,du1,du2,dv0,dv1,dv2;
% %   float D[3];
% %   float isect1[2], isect2[2];
% %   float du0du1,du0du2,dv0dv1,dv0dv2;
% %   short index;
% %   float vp0,vp1,vp2;
% %   float up0,up1,up2;
% %   float b,c,max;
% %
% %   /* compute plane equation of triangle(V0,V1,V2) */
% %   SUB(E1,V1,V0);
% %   SUB(E2,V2,V0);
% %   CROSS(N1,E1,E2);
% %   d1=-DOT(N1,V0);
% %   /* plane equation 1: N1.X+d1=0 */
% %
% %   /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
% %   du0=DOT(N1,U0)+d1;
% %   du1=DOT(N1,U1)+d1;
% %   du2=DOT(N1,U2)+d1;
% %
% %   /* coplanarity robustness check */
% % #if USE_EPSILON_TEST==TRUE
% %   if(fabs(du0)<EPSILON) du0=0.0;
% %   if(fabs(du1)<EPSILON) du1=0.0;
% %   if(fabs(du2)<EPSILON) du2=0.0;
% % #endif
% %   du0du1=du0*du1;
% %   du0du2=du0*du2;
% %
% %   if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute plane of triangle (U0,U1,U2) */
% %   SUB(E1,U1,U0);
% %   SUB(E2,U2,U0);
% %   CROSS(N2,E1,E2);
% %   d2=-DOT(N2,U0);
% %   /* plane equation 2: N2.X+d2=0 */
% %
% %   /* put V0,V1,V2 into plane equation 2 */
% %   dv0=DOT(N2,V0)+d2;
% %   dv1=DOT(N2,V1)+d2;
% %   dv2=DOT(N2,V2)+d2;
% %
% % #if USE_EPSILON_TEST==TRUE
% %   if(fabs(dv0)<EPSILON) dv0=0.0;
% %   if(fabs(dv1)<EPSILON) dv1=0.0;
% %   if(fabs(dv2)<EPSILON) dv2=0.0;
% % #endif
% %
% %   dv0dv1=dv0*dv1;
% %   dv0dv2=dv0*dv2;
% %
% %   if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute direction of intersection line */
% %   CROSS(D,N1,N2);
% %
% %   /* compute and index to the largest component of D */
% %   max=fabs(D[0]);
% %   index=0;
% %   b=fabs(D[1]);
% %   c=fabs(D[2]);
% %   if(b>max) max=b,index=1;
% %   if(c>max) max=c,index=2;
% %
% %   /* this is the simplified projection onto L*/
% %   vp0=V0[index];
% %   vp1=V1[index];
% %   vp2=V2[index];
% %
% %   up0=U0[index];
% %   up1=U1[index];
% %   up2=U2[index];
% %
% %   /* compute interval for triangle 1 */
% %   COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);
% %
% %   /* compute interval for triangle 2 */
% %   COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);
% %
% %   SORT(isect1[0],isect1[1]);
% %   SORT(isect2[0],isect2[1]);
% %
% %   if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
% %   return 1;
% % }
% %
% %
% % #define NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
% % { \
% %         if(D0D1>0.0f) \
% %         { \
% %                 /* here we know that D0D2<=0.0 */ \
% %             /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
% %                 A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
% %         } \
% %         else if(D0D2>0.0f)\
% %         { \
% %                 /* here we know that d0d1<=0.0 */ \
% %             A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
% %         } \
% %         else if(D1*D2>0.0f || D0!=0.0f) \
% %         { \
% %                 /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
% %                 A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
% %         } \
% %         else if(D1!=0.0f) \
% %         { \
% %                 A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
% %         } \
% %         else if(D2!=0.0f) \
% %         { \
% %                 A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
% %         } \
% %         else \
% %         { \
% %                 /* triangles are coplanar */ \
% %                 return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
% %         } \
% % }
% %
% %
% %
% % int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
% %                      float U0[3],float U1[3],float U2[3])
% % {
% %   float E1[3],E2[3];
% %   float N1[3],N2[3],d1,d2;
% %   float du0,du1,du2,dv0,dv1,dv2;
% %   float D[3];
% %   float isect1[2], isect2[2];
% %   float du0du1,du0du2,dv0dv1,dv0dv2;
% %   short index;
% %   float vp0,vp1,vp2;
% %   float up0,up1,up2;
% %   float bb,cc,max;
% %   float a,b,c,x0,x1;
% %   float d,e,f,y0,y1;
% %   float xx,yy,xxyy,tmp;
% %
% %   /* compute plane equation of triangle(V0,V1,V2) */
% %   SUB(E1,V1,V0);
% %   SUB(E2,V2,V0);
% %   CROSS(N1,E1,E2);
% %   d1=-DOT(N1,V0);
% %   /* plane equation 1: N1.X+d1=0 */
% %
% %   /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
% %   du0=DOT(N1,U0)+d1;
% %   du1=DOT(N1,U1)+d1;
% %   du2=DOT(N1,U2)+d1;
% %
% %   /* coplanarity robustness check */
% % #if USE_EPSILON_TEST==TRUE
% %   if(FABS(du0)<EPSILON) du0=0.0;
% %   if(FABS(du1)<EPSILON) du1=0.0;
% %   if(FABS(du2)<EPSILON) du2=0.0;
% % #endif
% %   du0du1=du0*du1;
% %   du0du2=du0*du2;
% %
% %   if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute plane of triangle (U0,U1,U2) */
% %   SUB(E1,U1,U0);
% %   SUB(E2,U2,U0);
% %   CROSS(N2,E1,E2);
% %   d2=-DOT(N2,U0);
% %   /* plane equation 2: N2.X+d2=0 */
% %
% %   /* put V0,V1,V2 into plane equation 2 */
% %   dv0=DOT(N2,V0)+d2;
% %   dv1=DOT(N2,V1)+d2;
% %   dv2=DOT(N2,V2)+d2;
% %
% % #if USE_EPSILON_TEST==TRUE
% %   if(FABS(dv0)<EPSILON) dv0=0.0;
% %   if(FABS(dv1)<EPSILON) dv1=0.0;
% %   if(FABS(dv2)<EPSILON) dv2=0.0;
% % #endif
% %
% %   dv0dv1=dv0*dv1;
% %   dv0dv2=dv0*dv2;
% %
% %   if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute direction of intersection line */
% %   CROSS(D,N1,N2);
% %
% %   /* compute and index to the largest component of D */
% %   max=(float)FABS(D[0]);
% %   index=0;
% %   bb=(float)FABS(D[1]);
% %   cc=(float)FABS(D[2]);
% %   if(bb>max) max=bb,index=1;
% %   if(cc>max) max=cc,index=2;
% %
% %   /* this is the simplified projection onto L*/
% %   vp0=V0[index];
% %   vp1=V1[index];
% %   vp2=V2[index];
% %
% %   up0=U0[index];
% %   up1=U1[index];
% %   up2=U2[index];
% %
% %   /* compute interval for triangle 1 */
% %   NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1);
% %
% %   /* compute interval for triangle 2 */
% %   NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);
% %
% %   xx=x0*x1;
% %   yy=y0*y1;
% %   xxyy=xx*yy;
% %
% %   tmp=a*xxyy;
% %   isect1[0]=tmp+b*x1*yy;
% %   isect1[1]=tmp+c*x0*yy;
% %
% %   tmp=d*xxyy;
% %   isect2[0]=tmp+e*xx*y1;
% %   isect2[1]=tmp+f*xx*y0;
% %
% %   SORT(isect1[0],isect1[1]);
% %   SORT(isect2[0],isect2[1]);
% %
% %   if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
% %   return 1;
% % }
% %
% % /* sort so that a<=b */
% % #define SORT2(a,b,smallest)       \
% %              if(a>b)       \
% %              {             \
% %                float c;    \
% %                c=a;        \
% %                a=b;        \
% %                b=c;        \
% %                smallest=1; \
% %              }             \
% %              else smallest=0;
% %
% %
% % inline void isect2(float VTX0[3],float VTX1[3],float VTX2[3],float VV0,float VV1,float VV2,
% % 	    float D0,float D1,float D2,float *isect0,float *isect1,float isectpoint0[3],float isectpoint1[3])
% % {
% %   float tmp=D0/(D0-D1);
% %   float diff[3];
% %   *isect0=VV0+(VV1-VV0)*tmp;
% %   SUB(diff,VTX1,VTX0);
% %   MULT(diff,diff,tmp);
% %   ADD(isectpoint0,diff,VTX0);
% %   tmp=D0/(D0-D2);
% %   *isect1=VV0+(VV2-VV0)*tmp;
% %   SUB(diff,VTX2,VTX0);
% %   MULT(diff,diff,tmp);
% %   ADD(isectpoint1,VTX0,diff);
% % }
% %
% %
% % #if 0
% % #define ISECT2(VTX0,VTX1,VTX2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1) \
% %               tmp=D0/(D0-D1);                    \
% %               isect0=VV0+(VV1-VV0)*tmp;          \
% % 	      SUB(diff,VTX1,VTX0);               \
% % 	      MULT(diff,diff,tmp);               \
% %               ADD(isectpoint0,diff,VTX0);        \
% %               tmp=D0/(D0-D2);
% % /*              isect1=VV0+(VV2-VV0)*tmp;          \ */
% % /*              SUB(diff,VTX2,VTX0);               \     */
% % /*              MULT(diff,diff,tmp);               \   */
% % /*              ADD(isectpoint1,VTX0,diff);           */
% % #endif
% %
% % inline int compute_intervals_isectline(float VERT0[3],float VERT1[3],float VERT2[3],
% % 				       float VV0,float VV1,float VV2,float D0,float D1,float D2,
% % 				       float D0D1,float D0D2,float *isect0,float *isect1,
% % 				       float isectpoint0[3],float isectpoint1[3])
% % {
% %   if(D0D1>0.0f)
% %   {
% %     /* here we know that D0D2<=0.0 */
% %     /* that is D0, D1 are on the same side, D2 on the other or on the plane */
% %     isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1);
% %   }
% %   else if(D0D2>0.0f)
% %     {
% %     /* here we know that d0d1<=0.0 */
% %     isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1);
% %   }
% %   else if(D1*D2>0.0f || D0!=0.0f)
% %   {
% %     /* here we know that d0d1<=0.0 or that D0!=0.0 */
% %     isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1);
% %   }
% %   else if(D1!=0.0f)
% %   {
% %     isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1);
% %   }
% %   else if(D2!=0.0f)
% %   {
% %     isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1);
% %   }
% %   else
% %   {
% %     /* triangles are coplanar */
% %     return 1;
% %   }
% %   return 0;
% % }
% %
% % #define COMPUTE_INTERVALS_ISECTLINE(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1,isectpoint0,isectpoint1) \
% %   if(D0D1>0.0f)                                         \
% %   {                                                     \
% %     /* here we know that D0D2<=0.0 */                   \
% %     /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
% %     isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
% %   }
% % #if 0
% %   else if(D0D2>0.0f)                                    \
% %   {                                                     \
% %     /* here we know that d0d1<=0.0 */                   \
% %     isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
% %   }                                                     \
% %   else if(D1*D2>0.0f || D0!=0.0f)                       \
% %   {                                                     \
% %     /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
% %     isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
% %   }                                                     \
% %   else if(D1!=0.0f)                                     \
% %   {                                                     \
% %     isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
% %   }                                                     \
% %   else if(D2!=0.0f)                                     \
% %   {                                                     \
% %     isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
% %   }                                                     \
% %   else                                                  \
% %   {                                                     \
% %     /* triangles are coplanar */                        \
% %     coplanar=1;                                         \
% %     return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
% %   }
% % #endif
% %
% % int tri_tri_intersect_with_isectline(float V0[3],float V1[3],float V2[3],
% % 				     float U0[3],float U1[3],float U2[3],int *coplanar,
% % 				     float isectpt1[3],float isectpt2[3])
% % {
% %   float E1[3],E2[3];
% %   float N1[3],N2[3],d1,d2;
% %   float du0,du1,du2,dv0,dv1,dv2;
% %   float D[3];
% %   float isect1[2], isect2[2];
% %   float isectpointA1[3],isectpointA2[3];
% %   float isectpointB1[3],isectpointB2[3];
% %   float du0du1,du0du2,dv0dv1,dv0dv2;
% %   short index;
% %   float vp0,vp1,vp2;
% %   float up0,up1,up2;
% %   float b,c,max;
% %   float tmp,diff[3];
% %   int smallest1,smallest2;
% %
% %   /* compute plane equation of triangle(V0,V1,V2) */
% %   SUB(E1,V1,V0);
% %   SUB(E2,V2,V0);
% %   CROSS(N1,E1,E2);
% %   d1=-DOT(N1,V0);
% %   /* plane equation 1: N1.X+d1=0 */
% %
% %   /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
% %   du0=DOT(N1,U0)+d1;
% %   du1=DOT(N1,U1)+d1;
% %   du2=DOT(N1,U2)+d1;
% %
% %   /* coplanarity robustness check */
% % #if USE_EPSILON_TEST==TRUE
% %   if(fabs(du0)<EPSILON) du0=0.0;
% %   if(fabs(du1)<EPSILON) du1=0.0;
% %   if(fabs(du2)<EPSILON) du2=0.0;
% % #endif
% %   du0du1=du0*du1;
% %   du0du2=du0*du2;
% %
% %   if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute plane of triangle (U0,U1,U2) */
% %   SUB(E1,U1,U0);
% %   SUB(E2,U2,U0);
% %   CROSS(N2,E1,E2);
% %   d2=-DOT(N2,U0);
% %   /* plane equation 2: N2.X+d2=0 */
% %
% %   /* put V0,V1,V2 into plane equation 2 */
% %   dv0=DOT(N2,V0)+d2;
% %   dv1=DOT(N2,V1)+d2;
% %   dv2=DOT(N2,V2)+d2;
% %
% % #if USE_EPSILON_TEST==TRUE
% %   if(fabs(dv0)<EPSILON) dv0=0.0;
% %   if(fabs(dv1)<EPSILON) dv1=0.0;
% %   if(fabs(dv2)<EPSILON) dv2=0.0;
% % #endif
% %
% %   dv0dv1=dv0*dv1;
% %   dv0dv2=dv0*dv2;
% %
% %   if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
% %     return 0;                    /* no intersection occurs */
% %
% %   /* compute direction of intersection line */
% %   CROSS(D,N1,N2);
% %
% %   /* compute and index to the largest component of D */
% %   max=fabs(D[0]);
% %   index=0;
% %   b=fabs(D[1]);
% %   c=fabs(D[2]);
% %   if(b>max) max=b,index=1;
% %   if(c>max) max=c,index=2;
% %
% %   /* this is the simplified projection onto L*/
% %   vp0=V0[index];
% %   vp1=V1[index];
% %   vp2=V2[index];
% %
% %   up0=U0[index];
% %   up1=U1[index];
% %   up2=U2[index];
% %
% %   /* compute interval for triangle 1 */
% %   *coplanar=compute_intervals_isectline(V0,V1,V2,vp0,vp1,vp2,dv0,dv1,dv2,
% % 				       dv0dv1,dv0dv2,&isect1[0],&isect1[1],isectpointA1,isectpointA2);
% %   if(*coplanar) return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
% %
% %
% %   /* compute interval for triangle 2 */
% %   compute_intervals_isectline(U0,U1,U2,up0,up1,up2,du0,du1,du2,
% % 			      du0du1,du0du2,&isect2[0],&isect2[1],isectpointB1,isectpointB2);
% %
% %   SORT2(isect1[0],isect1[1],smallest1);
% %   SORT2(isect2[0],isect2[1],smallest2);
% %
% %   if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
% %
% %   /* at this point, we know that the triangles intersect */
% %
% %   if(isect2[0]<isect1[0])
% %   {
% %     if(smallest1==0) { SET(isectpt1,isectpointA1); }
% %     else { SET(isectpt1,isectpointA2); }
% %
% %     if(isect2[1]<isect1[1])
% %     {
% %       if(smallest2==0) { SET(isectpt2,isectpointB2); }
% %       else { SET(isectpt2,isectpointB1); }
% %     }
% %     else
% %     {
% %       if(smallest1==0) { SET(isectpt2,isectpointA2); }
% %       else { SET(isectpt2,isectpointA1); }
% %     }
% %   }
% %   else
% %   {
% %     if(smallest2==0) { SET(isectpt1,isectpointB1); }
% %     else { SET(isectpt1,isectpointB2); }
% %
% %     if(isect2[1]>isect1[1])
% %     {
% %       if(smallest1==0) { SET(isectpt2,isectpointA2); }
% %       else { SET(isectpt2,isectpointA1); }
% %     }
% %     else
% %     {
% %       if(smallest2==0) { SET(isectpt2,isectpointB2); }
% %       else { SET(isectpt2,isectpointB1); }
% %     }
% %   }
% %   return 1;
% % }
% %
% %
% %
% %
% %








