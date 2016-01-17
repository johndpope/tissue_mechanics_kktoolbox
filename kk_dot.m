function [rvec] = kk_dot(u,v)
% u and v are assumed to be 3-vector(s). The columns being 3.
rvec = u(:,1).*v(:,1) + u(:,2).*v(:,2) + u(:,3).*v(:,3);