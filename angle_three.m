function A = angle_three(A,B,C)
% where A, B and C are 3-vectors: dimensions nx3
% returns the angle at A
u = B-A;
v = C-A;
A = acos(dot(u,v, 2));
