function [rvec] = kk_cross(u,v)
% u and v are assumed to be 3-vector(s). The columns being 3.
% % u = < a , b , c > and v = < d , e , f >
% % The cross product, noted by x, of the two vectors u and v given above is another vector w given by
% % w = u x v = < a , b , c > x < d , e , f > = < x , y , z >
% % with the components x, y and z given by:
% % x = b*f - c*e , y = c*d - a*f and z = a*e - b*d 

x = u(:,2).*v(:,3)-u(:,3).*v(:,2);
y = u(:,3).*v(:,1)-u(:,1).*v(:,3);
z = u(:,1).*v(:,2)-u(:,2).*v(:,1);
rvec = [x y z];