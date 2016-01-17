function [Xout Fout Vol] = rbf_remesh(dsites,F)
% based on radial basis functions interpolation
% uses iso2mesh library (alternatively matlab's poor isosurface)

%% configure
ep = 1; % Parameter for basis function
npu = 8;% Parameter for npu-by-npu-by-npu grid of PU cells
neval = 128;% Parameter for npu-by-npu-by-npu grid of PU cells
tol = 0.02;    % percentage tolerance on area and volume (for finding the correct object)
tol = 1.0;
opt.radbound = 2; % iso2mesh parameter
opt.distbound = 2; % iso2mesh parameter (see below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsites = double(dsites);
% calculate area, volume and normals of initial shape
disp('Calculating normals');
[A1, V1, v1, F_areas, h, H, wb, da, normals, K, k_g, dA] = triangulated_props(dsites, F, 0);
dfig(1);clf;plot_mesh(dsites,F);view(3);drawnow;
disp('Done');
bmino = min(dsites,[],1);  bmaxo = max(dsites,[],1);
bdimo = (bmaxo-bmino);
%%
fac = 1;
% Weight function for global Shepard partition of unity weighting
wf = @(e,r) r.^4.*(5*spones(r)-4*r);
% The RBF basis function for local RBF interpolation
rbf = @(e,r) r.^4.*(5*spones(r)-4*r);

% Load data points and compute bounding box
N = size(dsites,1);
bmin = min(dsites,[],1);  bmax = max(dsites,[],1);
bdim = max(bmax-bmin);
wep = npu/bdim;
% Add auxiliary points along normals "inside" and "outside"
% Find points with nonzero normal vectors and count them
withnormals = find(normals(:,1)|normals(:,2)|normals(:,3));
addpoints = length(withnormals);
% Distance along normal at which to place new points
delta = bdim/100;
% Create new points
dsites(N+1:N+addpoints,:) = ...
    dsites(withnormals,:) + delta*normals(withnormals,:);
dsites(N+addpoints+1:N+2*addpoints,:) = ...
    dsites(withnormals,:) - delta*normals(withnormals,:);
% Interpolant is implicit surface, i.e.,
% "original" points have rhs=0, "inside" rhs=-1, "outside" rhs=1
rhs = [zeros(N,1); fac *ones(addpoints,1); -fac * ones(addpoints,1)];
% Compute new bounding box
bmin = min(dsites,[],1);  bmax = max(dsites,[],1);
ctrs = dsites;
% Create neval-by-neval-by-neval equally spaced evaluation
% locations in bounding box
xgrid = linspace(bmin(1),bmax(1),neval);
ygrid = linspace(bmin(2),bmax(2),neval);
zgrid = linspace(bmin(3),bmax(3),neval);
[xe,ye,ze] = meshgrid(xgrid,ygrid,zgrid);
epoints = [xe(:) ye(:) ze(:)];
% Create npu-by-npu-by-npu equally spaced centers of PU cells
% in bounding box
puxgrid = linspace(bmin(1),bmax(1),npu);
puygrid = linspace(bmin(2),bmax(2),npu);
puzgrid = linspace(bmin(3),bmax(3),npu);
[xpu,ypu,zpu] = meshgrid(puxgrid,puygrid,puzgrid);
cellctrs = [xpu(:) ypu(:) zpu(:)];
cellradius = 1/wep;
% Compute Shepard evaluation matrix
DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);
SEM = wf(wep,DM_eval);
SEM = spdiags(1./(SEM*ones(npu^3,1)),0,neval^3,neval^3)*SEM;
% Build k-D trees for data sites and evaluation points
[tmp,tmp,datatree] = kdtree(dsites,[]);
[tmp,tmp,evaltree] = kdtree(epoints,[]);
Pf = zeros(neval^3,1); % initialize
for j=1:npu^3
    % Find data sites in cell j
    [pts,dist,idx] = kdrangequery(datatree,...
        cellctrs(j,:),cellradius);
    if (length(idx) > 0)
        % Build local interpolation matrix for cell j
        DM_data = DistanceMatrixCSRBF(dsites(idx,:),...
            ctrs(idx,:),ep);
        IM = rbf(ep,DM_data);
        % Find evaluation points in cell j
        [epts,edist,eidx] = kdrangequery(evaltree,...
            cellctrs(j,:),cellradius);
        % Compute local evaluation matrix
        DM_eval = DistanceMatrixCSRBF(epoints(eidx,:),...
            ctrs(idx,:),ep);
        EM = rbf(ep,DM_eval);
        % Compute local RBF interpolant
        localfit = EM * (IM\rhs(idx));
        % Accumulate global fit
        Pf(eidx) = Pf(eidx) + localfit.*SEM(eidx,j);
    end
end


%% we are done with generating the volume. Now we need to isosurface it.
d = 5;      % padding for a cleaner isosurfaced mesh
xgrid = linspace(bmin(1),bmax(1),neval + 2*d);
ygrid = linspace(bmin(2),bmax(2),neval + 2*d);
zgrid = linspace(bmin(3),bmax(3),neval + 2*d);
[xe,ye,ze] = meshgrid(xgrid,ygrid,zgrid);
Vol = padarray(reshape(Pf,neval,neval,neval),[d d d]);

%% uncomment for using matlab's poor implementation of isosurfacing (low quality mesh possible)
% % [F X] = isosurface(Vol,0);
% % facecell=finddisconnsurf(F);
%% uncomment for using vectorized marching cubes by Peter Hammer (low quality mesh possible)
% % [F,X,col] = MarchingCubes(xe,ye,ze,Vol,0);
% % facecell=finddisconnsurf(F);
%% uncomment for using iso2mesh library (low quality mesh possible, usually better than marching cubes)
%if method is 'cgalsurf' or 'cgalpoly':
%	     opt=a float number>1: max radius of the Delaunay sphere(element size)
%	     opt.radbound: same as above, max radius of the Delaunay sphere
%	     opt.distbound: maximum deviation from the specified isosurfaces
%	     opt(1,2,...).radbound: same as above, for each levelset
% opt.radbound = 2;
% opt.distbound = 2;
% % [X,F]=v2s(Vol,0, opt,'cgalsurf');
[X,elem, F]=v2m(Vol,0, 3,100);
 facecell=finddisconnsurf(F(:,1:3));

%% quickly determine which one of the objects is the correct one
indx = 0;
for(ix = 1:numel(facecell))
    [Xtemp,Ftemp]=removeisolatednode(X,facecell{ix});
    [Xtemp Ftemp] = fix_normals(X,F);
    bmin = min(Xtemp,[],1);  bmax = max(Xtemp,[],1);
    bdim_new = (bmax-bmin);
    Xtemp(:,1) = Xtemp(:,1)/bdim_new(1) * bdimo(1);
    Xtemp(:,2) = Xtemp(:,2)/bdim_new(2) * bdimo(2);
    Xtemp(:,3) = Xtemp(:,3)/bdim_new(3) * bdimo(3);
    [A, V, v2] = triangulated_props(Xtemp, Ftemp, 0);
    if(abs(A1-A)<(tol*A1) && abs(V1-V)<(tol*V1))
        Xout = Xtemp;
        Fout = Ftemp;%% now store the object that is relevant
        indx = ix;
        break;
    end
end
if(indx==0),
    warning('No legal object was generated');
    Xout = [];
    Fout = [];
end



% % dfig(ix + 1);clf;plot_mesh(Xtemp,Ftemp); view(3);drawnow;
% % disp('Comparison:');
% % disp([A1 A]);
% % disp([V1 V]);
% % disp([v1 v2]);

