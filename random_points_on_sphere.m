function [x y z t p] = random_points_on_sphere(n)
% generate a set of n*4 random points on the sphere
% method of Marsalgia 1972.

%% first quarter
 u0 = rand(2*n,1);
 u1 = rand(2*n,1);
 ind = (u0.^2 + u1.^2) <1;
 u0 = u0(ind);if length(u0)>n, u0 = u0(1:n);end
 u1 = u1(ind);if length(u1)>n, u1 = u1(1:n);end
 
 x1 = 2*u0.*sqrt(1-u0.^2-u1.^2);
 y1 = 2*u1.*sqrt(1-u0.^2-u1.^2);
 z1 = 1-2*(u0.^2+u1.^2);
 %% second quarter
 u0 = rand(2*n,1);
 u1 = rand(2*n,1);
 ind = (u0.^2 + u1.^2) <1;
 u0 = u0(ind);if length(u0)>n, u0 = u0(1:n);end
 u1 = u1(ind);if length(u1)>n, u1 = u1(1:n);end
 
 x2 = -2*u0.*sqrt(1-u0.^2-u1.^2);
 y2 = 2*u1.*sqrt(1-u0.^2-u1.^2);
 z2 = 1-2*(u0.^2+u1.^2);
  %% third quarter
 u0 = rand(2*n,1);
 u1 = rand(2*n,1);
 ind = (u0.^2 + u1.^2) <1;
 u0 = u0(ind);if length(u0)>n, u0 = u0(1:n);end
 u1 = u1(ind);if length(u1)>n, u1 = u1(1:n);end
 
 x3 = 2*u0.*sqrt(1-u0.^2-u1.^2);
 y3 = -2*u1.*sqrt(1-u0.^2-u1.^2);
 z3 = 1-2*(u0.^2+u1.^2);
   %% fourth quarter
 u0 = rand(2*n,1);
 u1 = rand(2*n,1);
 ind = (u0.^2 + u1.^2) <1;
 u0 = u0(ind);if length(u0)>n, u0 = u0(1:n);end
 u1 = u1(ind);if length(u1)>n, u1 = u1(1:n);end
 
 x4 = -2*u0.*sqrt(1-u0.^2-u1.^2);
 y4 = -2*u1.*sqrt(1-u0.^2-u1.^2);
 z4 = 1-2*(u0.^2+u1.^2);
 %%
 x = [x1(:);x2(:);x3(:);x4(:)];
 y = [y1(:);y2(:);y3(:);y4(:)];
 z = [z1(:);z2(:);z3(:);z4(:)];
 [t,p,r] = kk_cart2sph(x, y, z);