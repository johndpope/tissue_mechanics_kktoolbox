function [A d] = tri_spherical_prop(V1, V2, V3)
%%% calculate the triangle properties for a single spherical triangle
%%% Input: Vi containts theta and phi values that position the triangle on
%%% the unit sphere
%%% output: A: area
%%%         d: total perimeter length (if required)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin >1, X =[V1(:)';V2(:)';V3(:)'];
else
    X = V1;
end

[x y z] = kk_sph2cart(X(:,1), X(:,2),1);
Xc = [x(:) y(:) z(:)];
d = [orthodromic_distance(X(1,1),X(1,2), X(2,1),X(2,2)) ...
    orthodromic_distance(X(1,1),X(1,2), X(3,1),X(3,2)) ...
    orthodromic_distance(X(2,1),X(2,2), X(3,1),X(3,2))];
a = d(1);
b = d(2);
c = d(3);
s = a + b + c;
E = 4 * atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)));
A = abs(real(pi*E/180));

% % ang1 = acosd(kk_dot(Xc(1,:), Xc(2,:)));
% % ang2 = acosd(kk_dot(Xc(1,:), Xc(3,:)));
% % ang3 = acosd(kk_dot(Xc(2,:), Xc(3,:)));
% % A = pi * (ang1 + ang2 + ang3 - 180)/180;
% % A = abs(A);
%%%%
% d = [orthodromic_distance(X(1,1),X(1,2), X(2,1),X(2,2)) ...
%     orthodromic_distance(X(1,1),X(1,2), X(3,1),X(3,2)) ...
%     orthodromic_distance(X(2,1),X(2,2), X(3,1),X(3,2))];

%%
function d = orthodromic_distance(t1,p1, t2, p2)
%% all radians, assiming r = 1;
%% based on : http://en.wikipedia.org/wiki/Great-circle_distance

nominator = sqrt( (cos(t2)*sin(p2-p1)).^2 + (cos(t1).*sin(t2)-sin(t1).*cos(t2).*cos(p2-p1)).^2 );
denominator = sin(t1).*sin(t2) + cos(t1).*cos(t2).*cos(p2-p1);
d = atan2(nominator,denominator);


