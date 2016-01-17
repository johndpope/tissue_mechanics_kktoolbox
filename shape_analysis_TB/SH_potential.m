function [t,p] = SH_potential(t,p,niter, E,F, az, el)
%%% obtain a uniform configuration by using the derivative
%%% of spherical harmonic expansion of the sampling density

%% [1] Generate a fine "Eulerian" mesh and calcuate the average geodesic
%% distance between its nodes
dim = 80;P = partsphere(dim^2);[theta phi r] = kk_cart2sph(P(1,:)',P(2,:)',P(3,:)');
[u v w] = kk_sph2cart(theta,phi,1); M = [u v w];nE = size(M,1);
[C] = convhulln(M, {'Qt'}); % get a proper triangulation
counter = 0;
for ix = 1:nE,   % loop over the vertices
    fmemb = ismember(C, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
    [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
    links = [];
    for ik = 1:length(ig),links = [links C(ig(ik),:)];end
    L{ix} = unique(links(links~=ix));   llinks = L{ix};
    for jx = 1:length(llinks), counter = counter + 1;EE(counter,:) = [ix llinks(jx)];end;
end;
EE = sort(EE,2);      % sort the 2-vector according to its second dimension
EE = unique(EE,'rows');% this should reduce the size of E by 50%.
d = u(EE(:,1)).*u(EE(:,2)) + v(EE(:,1)).*v(EE(:,2)) + w(EE(:,1)).*w(EE(:,2));
dx = mean(acos(d));clear EE L C d P
%%%%%%% Prepare initial quantities
L_max = 6;
warning off MATLAB:divideByZero;
[L, K ] = indices_gen(1:(L_max + 1)^2);
M = length(L); %% number of functions in expansion
N = length(theta); %% number of data points
A  = zeros(N, M);
for S = 1:length(L),A(:, S) = ylk_cos_sin(L(S),K(S),phi', theta')';end

global P_LK
[Y_LK P_LK]		= precalc_ylk_cos_sin(p', t', L_max);
Y_LK_phi 		= (precalc_ylk_cos_sin_dphi(p', t', L_max));
[Y_LK_theta] 	= (precalc_ylk_cos_sin_dtheta(p', t', L_max));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save data_sh_mesh;


load data_sh_mesh;
nE = length(u);
step = .1;
temperature = 0.5;alph = 1;% dx = dx *10;
P = zeros(1,length(M)); % for the density potential
Def = zeros(1,length(M)); % for the deformation
Rt_def = 0;Rp_def = 0;
nV = length(t); [x y z] = kk_sph2cart(t,p,1); X = [x y z];
for iter = 1:20
    disp(iter);    
    %%%%  [pre 2] Calculate the deformation measure attached to each vertex
% %     def = [];
% %     for jx = 1:nV,  %loop over the vertices
% %         Ejx   = E(E(:,1)==jx,:);    % get only the E that we care about
% %         gdist = x(Ejx(:,1)).*x(Ejx(:,2)) + y(Ejx(:,1)).*y(Ejx(:,2)) + z(Ejx(:,1)).*z(Ejx(:,2));gdist = acos(gdist);
% %         def(jx) = sum((abs(gdist-dx)/dx));% measure of deformation
% %     end
    %%% [2] Spreading: assign each mesh point in M a density
    def_vec = [];
    for ix = 1:nE,  %% loop over the Eulerian grid
        P(ix) = 0;def_vec = [];
        for jx = 1:nV,  %loop over the vertices
            % calculate geodesic distance between ith M and jth P points
            r = u(ix).*x(jx) + v(ix).*y(jx) + w(ix).*z(jx);r = acos(r);
            if abs(r)<(5*dx),  
                P(ix) = P(ix)+ 1/4/dx*(1+cos(pi*r/2/dx));
% %                 def_vec = [def_vec def(jx)]; 
            end
        end
% %         if ~(isempty(def_vec)), Def(ix) = sum(def_vec)/length(def_vec);end
    end
    
    %%%%%%% [3] Expand the sampling density in SH
    %%%%%%% and [4] Calculate the first derivative of the sampling-density expansion
    P = P/norm(P)+temperature; [xD yD zD] = kk_sph2cart(theta',phi',P);
    R_pot = sqrt(xD.^2+yD.^2+zD.^2); clks = inv(A'*A)*(A'*R_pot');   %% Solve Linear Least squares
    clear cc ;cc(1,1,:) = clks;clks = cc(1,ones(1,length(p)),:);
    Rp_pot =(sum(clks.*Y_LK_phi,3));Rt_pot = (sum(clks.*Y_LK_theta,3));
    Rt_pot(1) = 0;Rt_pot(end) = 0;
    %%%%%%% [3] Expand the local deformation
    %%%%%%% and [4] Calculate the first derivative of the local deformation
% %     Def = Def-norm(Def)+1; [xD yD zD] = kk_sph2cart(theta',phi',Def);
% %     [clks] = sh_analysis_LS([xD yD zD], L_max, 1);
% %     
% %     R_def = sqrt(xD.^2+yD.^2+zD.^2); clks = inv(A'*A)*(A'*R_def');   %% Solve Linear Least squares
% %     clear cc ;cc(1,1,:) = clks;clks = cc(1,ones(1,length(p)),:);
% %     Rp_def =(sum(clks.*Y_LK_phi,3));Rt_def = (sum(clks.*Y_LK_theta,3));
% %     Rt_def(1) = 0;Rt_def(end) = 0;
    
    %%%%% [5] Take the step
    t = t - (Rt_pot' + alph*Rt_def') * step; 
    p = p - (Rp_pot' + alph*Rp_def') * step;
    [x y z] = kk_sph2cart(t,p,1); X = [x y z];
    %%%%% [6] Plot
    plot_state(x,y,z,F, iter, 1,az, el);
end
save tp t p