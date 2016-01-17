function [c, ceq] = constraints(X,flag)
% global variables
global F nvert F_areas_relative pcount DM invDM
t = X(1:(nvert)); p = X(nvert+1:end); p = mod(p,2*pi);[u, v, w] = kk_sph2cart(t,p ,1);% convert coordinates (on the sphere) to the Cartesian
x1 = u(F(:,1)); y1 = v(F(:,1));z1 =  w(F(:,1));x2 = u(F(:,2)); y2 = v(F(:,2));z2 =  w(F(:,2));x3 = u(F(:,3)); y3 = v(F(:,3));z3 =  w(F(:,3));
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2); twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
F_areas = twoA./2;
%%% inequality constraints
nq = sqrt(sum(q.*q,2));nr = sqrt(sum(r.*r,2));
ttr = acos(dot(q,r,2)./nq./nr)*180/pi;
fac = 5; 
c = [ttr-(180-fac); -ttr-fac];
% c = []
%%%% equality constraints
ceq= [];