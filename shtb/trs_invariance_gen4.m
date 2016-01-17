function [X_o] = trs_invariance_gen4(X_o)
%%% Calculate the Euler angles a,b,g needed for a L = 1 based rotational invariant
%%% shape description of X_o shape.

[xc yc zc] = get_xyz_clks(times_Basis(X_o));
% translation invariance
xc(1) = 0;yc(1) = 0;zc(1) = 0;

Cr = [xc(:) yc(:) zc(:)]';


% % parameter space rotation
Ax = xc(2:4);Ay = yc(2:4);Az = zc(2:4);
A = [Ax(:) Ay(:) Az(:)];
[Ruvw D] = eig(A'*A);

Cr = Ruvw'*Cr;
xc = Cr(1,:);yc = Cr(2,:);zc = Cr(3,:); 

% object space rotation
A = [xc(2:4)';yc(2:4)';zc(2:4)'];
%%% diagonalize the ellipsoid
[V D] = eig(A'*A);
Cr = V'*Cr;  % by applying the transformation to the L = 1 coefficients (and don't forget weighting by the basis normalization) we obtain the diagonalized form
xc = Cr(1,:);yc = Cr(2,:);zc = Cr(3,:); 
X_o = over_Basis([xc(:)' yc(:)' zc(:)']);
