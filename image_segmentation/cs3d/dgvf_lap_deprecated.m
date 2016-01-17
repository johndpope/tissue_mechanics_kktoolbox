function [un, vn, wn] = dgvf_lap(I, niter, miu, dt, dx, dy, dz)
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

verbose = 1;
[un, vn, wn] = gradient(double(I)); 
E = -sqrt(un.^2 + vn.^2 + wn.^2); 		% Energy as in Xu and Prince 1998 Equation 2
                                        % if the image is a line drawing
                                        % then use: E = I --> Eq. 4
[fx, fy, fz] = gradient(-E); clear E;	
b = fx.^2 + fy.^2 + fz.^2;              % Equation 15
c{1} = fx; clear fx
c{2} = fy; clear fy
c{3} = fz; clear fz

r = miu*dt/dx/dy/dz;                	% Equation 17 (Xu and Prince 1998)
if (dt>(dx*dy*dz/6/miu)), disp('convergence not guaranteed!!!! Resetting dt');dt =(dx*dy*dz/6/miu/2);disp(dt);end; % Equation 18

vec_n{1} = un;
vec_n{2} = vn;
vec_n{3} = wn;

if verbose, disp('Calculating 3D diffusion gradient vector field using prallel toolbox');end
parfor dim = 1:3       % parallelize over three processors if possible
        vec_n{dim} = time_steps(vec_n{dim}, miu, b, c{dim}, niter);
end
un = vec_n{1};
vn = vec_n{2};
wn = vec_n{3};

if sum(isnan([un(:)' vn(:)' wn(:)']));error('diffusion gradient calculation failed');end

%%
function un = time_steps(f, miu, b, c, niter)
%% Start the time stepping (Xu Prince 1998 p.363 using 3D version of Equation 16 (Xu and Prince 1998)
for n = 1:niter,  
  [N M O] = size(f);
 
    xi = 2:M-1;
    yi = 2:N-1;
    zi = 2:O-1;
     % Coners
    f([1 N], [1 M], [1 O]) = f([3 N-2], [3 M-2], [3 O-2]);
     % Edges
    f([1 N], [1 M], zi) = f([3 N-2], [3 M-2], zi);
    f(yi, [1 M], [1 O]) = f(yi, [3 M-2], [3 O-2]);
    f([1 N], xi, [1 O]) = f([3 N-2], xi, [3 O-2]);
     % Faces
    f([1 N], xi, zi) = f([3 N-2], xi, zi);
    f(yi, [1 M], zi) = f(yi, [3 M-2], zi);
    f(yi, xi, [1 O]) = f(yi, xi, [3 O-2]);   
    
    f = f + miu*6*del2(f)-(f-c).*b;
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
