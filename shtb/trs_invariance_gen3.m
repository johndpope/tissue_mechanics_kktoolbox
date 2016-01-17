function [X_o] = trs_invariance_gen3(X_o)
%%% Calculate the Euler angles a,b,g needed for a L = 1 based rotational invariant
%%% shape description of X_o shape.
X_o = nocs2cs(X_o);
[xc yc zc] = get_xyz_clks(X_o);

Ax = xc(2:4);[Vx Dx] = eig(Ax'*Ax);
Ay = yc(2:4);[Vy Dy] = eig(Ay'*Ay);
Az = zc(2:4);[Vz Dz] = eig(Az'*Az);
V = Vx*Vy*Vz;



M = [Ax(:) Ay(:) Az(:)];
[V D] = eig(M'*M);


A = [xc(:) yc(:) zc(:)]';
Ar = V'*A;
xc = Ar(1,:);yc = Ar(2,:);zc = Ar(3,:);







X_o = cs2nocs([xc(:)' yc(:)' zc(:)']);
