function [ai, bi, dAi] = triangle_strain(T1, T2)
% returns the strains when deforming triangle T1 to become T2
% Input:
%       T1 and T2 are each 3-vectors with 3 rows, with columns corresponding to
%       x y and z coordinates. Each row is a vertex in the triangle in the
%       sequence V1, V2 and V3
% The algorithm is based on SoftMatter article Khairy et al 2012, Supp note 7 and Figure S-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1 = [];
l2 = [];
V1 = T1(1,:);
V2 = T1(2,:);
V3 = T1(3,:);
v1 = T2(1,:);
v2 = T2(2,:);
v3 = T2(3,:);
dm12 = V2-V1;
dm13 = V3-V1;
ds12 = v2-v1;
ds13 = v3-v1;

%%%% according to Lim Wortis 2006 
lxly = norm(cross(v1, v2))/norm(cross(V1, V2));
lxsq = (norm(v1))^2/(norm(V1))^2;
P = (dot(V1, V2)/(norm(V1))^2 - dot(v1, v2)/(norm(v1))^2)^2 * (norm(v1))^2;
Q = (norm(v2 - (dot(v1,v2)/(norm(v1))^2)*v1))^2;
tanphisq = P/Q;
%%%%
ai = lxly-1;
bi = 1/2 * (lxsq/lxly + lxly/lxsq * (1 + tanphisq)-2);
dAi = 1/2 * (norm(cross(V1, V2)));


