function [xclks,yclks, zclks, A, V, v, da,CI] = volume_gen(X,xclks, yclks, zclks,x_pixel, y_pixel,z_pixel,t,p,Y_LK,F, PSF, RAW, L_max_fit, plotflag)
%% Optimize the shape parameters using the convolved volume%%
global Y_LK PSF_fft xshift yshift zshift x_pixel y_pixel z_pixel RAW F normXo plotflag 
global fac frame_vec Xmx Xmy Xmz Iobj Iobj2 chisq itercount Xo ind_RAW ind_parm full_parm
global Y_LK_g P_LK Y_LK_phi Y_LK_theta P_LK_T Y_PP Y_TT Y_TP wp wt st ct stsq sp spsq cp cpsq ctsq ctcu ctqad r_o
global R Y_LKc C wbvec_x wbvec_y wbvec_z


xclks_start = xclks;yclks_start = yclks;zclks_start = zclks;
clc;xdim = size(PSF, 1);ydim = size(PSF, 2);zdim = size(PSF, 4);
PSF_fft = fftn(squeeze(PSF), [xdim ydim zdim]);clear PSF;
RAW = mat2gray(RAW, [min(min(RAW(:,:,1,32))) max(max(RAW(:,:,1,32)))]);
for ix = 1:size(RAW,4),RAW(:,:,1,ix) = medfilt2(RAW(:,:,1,ix),[3 3]);end;
montage(mat2gray(RAW));title('raw data');drawnow
thresh = 0.07;ind_RAW = find(RAW(:)>thresh);ind_RAW = find(RAW(:));

full_parm = [xclks(:);yclks(:);zclks(:)];
nparms = (L_max_fit+1)^2;nc = length(xclks);
ind_parm = [1:nparms nc+1:nc+nparms 2*nc+1:2*nc+nparms];
Xo = full_parm(ind_parm); %full_parm(full_parm<) = 0;disp(sum(logical(full_parm)));
normXo = 1;%norm(Xo);  %  will be used for scaling
Xo = Xo/normXo;
disp(['Using ' num2str(length(Xo)) ' parameters total']);figure
plotflag = 1;
%% prepare initial quantities for chisquare calculation
%fac = 1/1;frame_vec = 17:3:57;data_vec = RAW(:,:,1,frame_vec);data_vec = data_vec(:);
Iobj = zeros(xdim, ydim, 1, zdim);Iobj2 = Iobj; chisq = RAW(ind_RAW)*0;itercount = 1;disp(['Length of residual vector: ' num2str(length(chisq))]);
Xmx = zeros(size(Y_LK, 1), 1);Xmy = Xmx;Xmz = Xmx;

%% prepare initial quantities for single mode bending energy calculation
gdim = 30;r_o = 20;L_max = round(sqrt(length(xclks))-1);
R = zeros(gdim*gdim,1);
[t_g wt]                = gaussquad(gdim, 0, pi);
[p_g wp]                = gaussquad(gdim,0,2*pi);
[p_g t_g]               = meshgrid(p_g,t_g);
[wp wt]                 = meshgrid(wp, wt);
[Y_LK_g P_LK]			= precalc_ylk_cos_sin(p_g, t_g, L_max);
Y_LK_phi 				= precalc_ylk_cos_sin_dphi(p_g, t_g, L_max);
[Y_LK_theta P_LK_T] 	= precalc_ylk_cos_sin_dtheta(p_g, t_g, L_max);
Y_PP 					= precalc_ylk_cos_sin_dphiphi(p_g, t_g, L_max);
Y_TT 					= precalc_ylk_cos_sin_dthetatheta(p_g, t_g, L_max);
Y_TP 					= precalc_ylk_cos_sin_dthetaphi(p_g, t_g, L_max);

wt = reshape(wt,gdim*gdim,1);wp = reshape(wp,gdim*gdim,1);
p_g = reshape(p_g,gdim*gdim,1);t_g = reshape(t_g,gdim*gdim,1);
Y_LK_g = reshape(Y_LK_g,gdim*gdim,(L_max+1)^2);
Y_LK_theta = reshape(Y_LK_theta,gdim*gdim,(L_max+1)^2);
Y_LK_phi = reshape(Y_LK_phi,gdim*gdim,(L_max+1)^2);
Y_TT = reshape(Y_TT,gdim*gdim,(L_max+1)^2);
Y_PP = reshape(Y_PP,gdim*gdim,(L_max+1)^2);
Y_TP = reshape(Y_TP,gdim*gdim,(L_max+1)^2);
st  = sin(t_g);sp = sin(p_g); ct = cos(t_g); cp = cos(p_g);
ctsq = ct.^2; cpsq = cp.^2; stsq = st.^2; spsq = sp.^2;
ctcu = ct.^3; ctqad = ct.^4;

%%%%%%%%%% Precalculate necessary quantities for the fast object curvature calculation
global C HidAi dAi HidAi_mx dAij V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix u v w crossqpr

[Areac,Vc,vc,tc,pc,Xc,C,Y_LKc]=plot_sh(xclks, yclks, zclks, 80);% Calculate shape properties based on triangulation and plot
[E,L] = edge_info(Xc,C);save data_intermediate1
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
X_o = [xclks; yclks; zclks]';
disp(' Bending energy calculation initialized');
% clear nedges phi theta tr tr2 wp wt 
% clear x y z xclks yclks zclks sing r llinks links jx jg ix ik
% clear ig h gdim g2 g1 fmemb ememb edg counter c V_far 
% clear V_e VT VE_tr1 VE_tr2 V LEV E A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make clks for a sphere and calculate the perturbation bending energy vectors
[x, y, z] = mesh_gen(gdim, 1);nvert = length(X);
[t p Radius] = kk_cart2sph(x, y, z);
[SA, SV, Sh, Sv, Sxclks, Syclks, Szclks] = sh_projection2(gdim, L_max, L_max, L_max, x, y, z, t, p);
[Bx, By, Bz, wbvec_x, wbvec_y, wbvec_z] =  single_mode_bending_energy(Sxclks(:),Syclks(:),Szclks(:));

save data_initialization_volume_gen;
%[B] =  bending_energy(xclks(:),yclks(:),zclks(:));                              % The smoothing functional
%[Fun_vec] = volume_gen_objective_vec(Xo);drawnow;


%%%%%
warning off;
% save data_volume_gen_intermediate;
load data_initialization_volume_gen;
Xo = full_parm(ind_parm);
disp('Optimization started');
tic
options = optimset('MaxFunEvals', 80000,'DiffMaxChange', 1.0,'DiffMinChange', 1e-2,'FunValCheck', 'on',...
            'DerivativeCheck','off','GradObj','off','GradConstr','off',...
            'MaxIter', 15,'Display', 'iter','Diagnostics', 'off','LevenbergMarquardt','on',...
            'TolCon', 1e-5,'TolFun', 1e-4,'TolX', 1e-4, 'MaxSQPIter',30,...
            'LargeScale','off', 'LineSearchType','quadcubic',...
            'PrecondBandWidth',5,'GoalsExactAchieve',0, 'MeritFunction','singleobj',...%'TypicalX',0 * ones(length(Xo)+1,1),...
            'MaxPCGIter','max(1,floor(numberOfVariables/2))');        
warning off;
%[X, resnorm,residual,exitflag,output,lambda, Jacobian] = lsqnonlin(@volume_gen_objective,Xo,[],[],options);warning on;
% lambdavec = 1*10.^(-(:20));
% for ix = 9:(9+length(lambdavec))
%      lambda = lambdavec(ix);
lambda = 10^(-2.2);
    Xo = [xclks(:); yclks(:); zclks(:)];
    % goal = [200 B lambda*Bx lambda*By lambda*Bz];weight = abs(goal);
    % [X, fval, attainfactor] = fgoalattain(@volume_gen_objective_vec,Xo,goal, weight,[],[],[],[],[],[],[],options);

    %[X,fval,exitflag,output] = fminsearch(@volume_gen_objective_amoeba,Xo,options,lambda);
    [X,fval,exitflag,output,grad] = fminunc(@volume_gen_objective_amoeba,Xo,options,lambda);
    toc

    % clks = X(1:length(X)-3);
    clks = X;
    nc = round((length(clks))/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
    %nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
    [fval,Eb,chisq] = volume_gen_objective_amoeba([xclks(:); yclks(:); zclks(:)], lambda);
    plot_sh(xclks,yclks,zclks);
    
    xclks_vec(:,ix) = xclks;
    yclks_vec(:,ix) = yclks;
    zclks_vec(:,ix) = zclks;
    fval_vec(ix) = fval;
    Eb_vec(ix) = Eb;
    chisq_vec(ix) = chisq;
    disp(lambda);
    disp([fval Eb chisq]);
% end

save data_volume_gen_after_fitting;
%load data_volume_gen_after_fitting
%%%% calculating 95% confidence intervals on the estimated parameters
%%%% Note that the length of residual must be larger than length of clks
% if length(X)<length(residual),
%     CI = nlparci(X, residual, Jacobian);
% end
% save data_volume_gen_after_fitting;

%%% Calculate shape properties of optimized shape
 load data_volume_gen_after_fitting;
% [R] = volume_gen_objective_vec(X);
A = 0;V = 0;v = 0;da = 0;

X = X.*normXo;clks = full_parm;clks(ind_parm) = X(1:end);
% CI = CI * normXo;
%xshift = X(end-2);yshift = X(end-1);zshift = X(end);
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
% xCI = CI(1:nc,:);yCI = CI(nc+1:2*nc,:);zCI = CI(2*nc+1:end-3,:);
% gdim = size(Y_LK,1);
% Xm = [sum(xclks(ones(gdim,1),:).*Y_LK,2)  sum(yclks(ones(gdim,1),:).*Y_LK,2) sum(zclks(ones(gdim,1),:).*Y_LK,2)];
% Xm = [Xm(:,1)+xshift Xm(:,2)+yshift Xm(:,3)+zshift];

Xmx(:) = Y_LK*xclks(:) + xshift;
Xmy(:) = Y_LK*yclks(:) + yshift;
Xmz(:) = Y_LK*zclks(:) + zshift;
Xm = [Xmx Xmy Xmz];
%%%%%%% Plot the results
close all;

% figure; errorbar(1:length(X), X, abs(X - CI(:,1)), abs(X - CI(:,2)));
%figure; plot(X(1:end-3),'b-');hold on;plot(CI(1:end-3,1),'r-');plot(CI(1:end-3,2),'r-');hold off;
figure;montage(mat2gray(Iobj));figure; montage(mat2gray(RAW));drawnow;

% figure;subplot(3,1,1);plot(xclks_start,'k');hold on;plot(xclks, 'b');hold on;plot(xCI(:,1),'r-');plot(xCI(:,2),'r-');hold off;
% subplot(3,1,2);plot(yclks_start,'k');hold on;plot(yclks, 'b');hold on;plot(yCI(:,1),'r-');plot(yCI(:,2),'r-');hold off;
% subplot(3,1,3);plot(zclks_start,'k');hold on;plot(zclks, 'b');hold on;plot(zCI(:,1),'r-');plot(zCI(:,2),'r-');hold off;

figure;subplot(2,1,1);plot_sh(xclks_start,yclks_start,zclks_start);subplot(2,1,2);plot_sh(xclks,yclks,zclks);
patch('Vertices', Xm, 'Faces', F,'FaceColor', 'red','EdgeColor','none');daspect([1 1 1]);axis off;
light; lighting gouraud;        view(3);        [A, V, v_red] = triangulated_props(Xm, F,0);
str = sprintf('Area: %.4f    Vol.: %.4f  v: %.4f   iter:%d',A,-V,-v_red,itercount);title(str);drawnow;





