function [un, vn, wn] = dgvf_calc(I, niter, miu, dt, dx, dy, dz, un_last, vn_last, wn_last)
% Generate the diffusion gradient vector field as in Xu and Prince 1998
% dgvf_calc is the three dimensional extension of the 2D version described in Equation 12
% Xu and Prince 1998,"Snakes, Shapes, and Gradient Vector Flow", IEEE Transactions on Image Processing Vol.7(3)
% Input:
%			I: three dimensional image (matrix)
%			miu: the smoothing parameter (more smoothing -- higher miu-- for noisy images)
%			niter: number of iterations
%			dt: time step
%			dx, dy,dz: pixel spacing, set all to one for isotropic images
% Output:
%			un, vn, wn: the gradient vector field
% Example usage:
%			[un vn wn] = dgvf_calc(I,sqrt(numel(I)), 0.5, 1, 1, 1, 1)
% Note: circshift has been pointed out to be slow. Alternatives are
% required.
% Author:   Dr.Khaled Khairy, Janelia Farm Research Campus, Howard Hughes Medical Institute.
%			June 2011. Please send corrections or comments to khairyk@janelia.hhmi.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % global I1m1 I11 I2m1 I21 I3m1 I31 real

%%% configure
verbose = 0;
GPU = 0;
PRLL = 1;
if nargin==7, PRLL = 1;end
use_matlab_gradient = 0;        % there is really no reason to use matlab's slower gradient method.


if verbose, disp('Calculating gradient of image.');end
if use_matlab_gradient
    [un, vn, wn] = gradient((I));
else
    un = DGradient(I, [], 2);
    vn = DGradient(I, [], 1);
    wn = DGradient(I, [], 3);
end
clear I
E = -sqrt(un.^2 + vn.^2 + wn.^2); 		% Energy as in Xu and Prince 1998 Equation 2% if the image is a line drawing then use: E = I --> Eq. 4

if verbose, disp('Calculating gradient of force.');end
if use_matlab_gradient
    [fx, fy, fz] = gradient(-E);
else
    fx = DGradient(-E, [], 2);
    fy = DGradient(-E, [], 1);
    fz = DGradient(-E, [], 3);
end
clear E;
    
b = fx.^2 + fy.^2 + fz.^2;              % Equation 15
c{1} = (b.*fx); clear fx
c{2} = (b.*fy); clear fy
c{3} = (b.*fz); clear fz

r = miu*dt/dx/dy/dz;                	% Equation 17 (Xu and Prince 1998)
if (dt>(dx*dy*dz/6/miu)), disp('convergence not guaranteed!!!! Resetting dt');dt =(dx*dy*dz/6/miu/2);disp(dt);end; % Equation 18
vec_n = cell(3,1);

if nargin == 10     % then use the gradient field given
    clear un vn wn;
    vec_n{1} = (un_last);
    vec_n{2} = (vn_last);
    vec_n{3} = (wn_last);
else                % just use the gradient of the image as initialization
    vec_n{1} = (un);
    vec_n{2} = (vn);
    vec_n{3} = (wn);
end
b = (b);
dt = (dt);
r = (r);
bdt = (b.*dt);


if GPU
    if verbose, disp('Calculating 3D diffusion gradient vector field using GPU');end
    gcache flush;
    dt = gsingle(dt);
    r  = gsingle(r);
    niter = gsingle(niter);
    gsync
    for dim = 1:3       
       if verbose, fprintf('dimension: %d\n', dim);end
        A = gsingle(vec_n{dim});
        B = gsingle(bdt);
        C = gsingle(c{dim}*dt);
        vec_n{dim} = time_steps(A, B, r, C , niter);
    end
    gsync
else
    if verbose, disp('Calculating 3D diffusion gradient vector field');end
    if matlabpool('size')>=3 && PRLL
        fprintf('\nUsing parallel toolbox to calculate dgvf\n');
        parfor dim = 1:3       % parallelize over three processors
            if verbose, fprintf('dimension: %d\n', dim);end
            vec_n{dim} = time_steps((vec_n{dim}), bdt, r, (c{dim})*dt, niter);
        end
    else
        fprintf('Using single CPUs to calculate dgvf serially');
        for dim = 1:3       % then just do it sequentially
            if verbose, fprintf('dimension: %d\n', dim);end
            vec_n{dim} = time_steps((vec_n{dim}), bdt, r, (c{dim})*dt, niter);
        end
    end
end

un = double(vec_n{1});
vn = double(vec_n{2});
wn = double(vec_n{3});
if sum(isnan([un(:)' vn(:)' wn(:)']));error('diffusion gradient calculation failed');end

%%
function un = time_steps(un, bdt, r, cdt, niter)
%% Start the time stepping (Xu Prince 1998 p.363 using 3D version of Equation 16 (Xu and Prince 1998)
%%%% speed up suggestion: replace circshift with a custom written function
%%%% that avoids the error checks. 

persistent un1
if isempty(un1)
    un1 = zeros(size(un));
end

thresh = 1e-10;
% cutoff = thresh*numel(un);
cutoff = 0.000209;
fprintf('\nUsing cutoff %.6f\n', cutoff);
%%% uncomment to use without pre-indexing
counter = 1;
convergence = 0;
convergence_test = 10;
while counter<niter && (convergence==0)
    un1 = (1-bdt).*un+ r.*(  ...
          circshift(un,[-1  0  0])...
        + circshift(un,[ 1  0  0])...
        + circshift(un,[ 0 -1  0] )...
        + circshift(un,[ 0  1  0] )...
        + circshift(un,[ 0  0 -1] )...
        + circshift(un,[ 0  0  1])...
        - 6.*un)...
        + cdt;
    
    %% % sosi----
%     vec = 1:7;
%     unxm = circshift(un,[-1  0  0]);
%     unxp = circshift(un,[1  0  0]);
%     unym = circshift(un,[0  -1  0]);
%     unyp = circshift(un,[0  1  0]);
%     unzm = circshift(un,[0  0  -1]);
%     unzp = circshift(un,[0  0  1]);
%     disp([un(vec)' unxm(vec)' unxp(vec)' unym(vec)' unyp(vec)' unzm(vec)' unzp(vec)']);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if mod(counter,convergence_test)==0,
%         res = sum((un1(:)-un(:)).^2);
%         %disp(res);
%         if res<cutoff, convergence = 1;end
%     end
    un(:) = un1(:);
    counter = counter + 1;
end

if counter<niter,
    disp('DGVF converged');
else
    disp('DGVF reached maximum number of iterations.');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Start the time stepping (Xu Prince 1998 p.363 using 3D version of Equation 16 (Xu and Prince 1998)
% for n = 1:niter,
%     if verbose,disp(['Diffusion gradient calculation iteration : ' num2str(n) '  of  ' num2str(niter)]);end
%
%     un = (1-b.*dt).*un+ r.*(  circshift(un,[-1  0  0])...
%         + circshift(un,[ 1  0  0])...
%         + circshift(un,[ 0 -1  0] )...
%         + circshift(un,[ 0  1  0] )...
%         + circshift(un,[ 0  0 -1] )...
%         + circshift(un,[ 0  0  1])...
%         - 6.*un)...
%         + c1.*dt;
%
%     vn = (1-b.*dt).*vn+ r.*(  circshift(vn,[-1  0  0])...
%         + circshift(vn,[ 1  0  0])...
%         + circshift(vn,[ 0 -1  0] )...
%         + circshift(vn,[ 0  1  0] )...
%         + circshift(vn,[ 0  0 -1] )...
%         + circshift(vn,[ 0  0  1])...
%         - 6.*vn)...
%         + c2.*dt;
%     wn = (1-b.*dt).*wn+ r.*(  circshift(wn,[-1  0  0])...
%         + circshift(wn,[ 1  0  0])...
%         + circshift(wn,[ 0 -1  0] )...
%         + circshift(wn,[ 0  1  0] )...
%         + circshift(wn,[ 0  0 -1] )...
%         + circshift(wn,[ 0  0  1])...
%         - 6.*wn)...
%         + c3.*dt;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
