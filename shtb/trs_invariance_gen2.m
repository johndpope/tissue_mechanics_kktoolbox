function [X_o] = trs_invariance_gen2(X_o)
%%% Calculate the Euler angles a,b,g needed for a L = 1 based rotational invariant
%%% shape description of X_o shape.
X_o = nocs2cs(X_o);
[xc yc zc] = get_xyz_clks(X_o);

Ax = xc(2:4);[Vx Dx] = eig(Ax'*Ax);
Ay = yc(2:4);[Vy Dy] = eig(Ay'*Ay);
Az = zc(2:4);[Vz Dz] = eig(Az'*Az);
V = Vx*Vy*Vz;
A = [xc(:) yc(:) zc(:)]';
Ar = V'*A;
xc = Ar(1,:);yc = Ar(2,:);zc = Ar(3,:);

% xc(2:4) = V'*xc(2:4)';
% yc(2:4) = V'*yc(2:4)';
% zc(2:4) = V'*zc(2:4)';
% 
% Ay = yc(2:4);
% [V D] = eig(Ay'*Ay);
% xc(2:4) = V'*xc(2:4)';
% yc(2:4) = V'*yc(2:4)';
% zc(2:4) = V'*zc(2:4)';
% 
% 
% Az = zc(2:4);
% [V D] = eig(Az'*Az);
% xc(2:4) = V'*xc(2:4)';
% yc(2:4) = V'*yc(2:4)';
% zc(2:4) = V'*zc(2:4)';

% a = 0.5;    % rotate around y
% b = 0.0;    % rotate around x
% g = 0.0;    % rotate around y
% xc = sh_rot(xc,g,b,a);
% yc = sh_rot(yc,g,b,a);
% zc = sh_rot(zc,g,b,a);

X_o = [xc(:)' yc(:)' zc(:)'];


% A = [xc(2:4)' yc(2:4)' zc(2:4)']';  % coefficient matrix of L = 1
% [V,D] = eig(A'*A,'nobalance');                  % calculate eigenvectors and eigenvalues
% Ap = A*V;                            % Rotation matrix V that diagonalizes A is applied to A
% C  = [xc(:) yc(:) zc(:)]';
% Cp = C'*V;



% % % %% Step 1: Maximize the trace of the L = 1 associated tensor matrix 
% % % g = 0.5;b = -0.5;a = 0.0; 
% % % ang_o = [g b a]; %starting values
% % % [X] = rotation_objective(ang_o,xc,yc,zc);
% % % options =   optimset('NonlEqnAlgorithm ', 'gn', 'LineSearchType', 'cubicpoly', ...
% % %             'MaxFunEvals', 5e3,'MaxIter', 30,'DiffMaxChange', .1,'DiffMinChange', 1e-3,...
% % %             'DerivativeCheck','off','LargeScale','off','GradObj','off','GradConstr','off',...
% % %             'Display', 'off','Diagnostics', 'off',...
% % %             'TolCon', 1e-2,'TolFun', 1e-8,'TolX', 1e-10,'MaxSQPIter',100);
% % %         
% % % warning off;[ang] = fsolve(@rotation_objective,ang_o,options,xc,yc,zc);warning on;
% % % g = ang(1);b = ang(2);a = ang(3);
% % % 
% % % %% Step 2: Rotate the coefficients (shape); rotate the clks to align with the major axis
% % % xrot = xc;yrot = yc;zrot = zc;
% % % [xrot] = sh_rot(xrot, g, b, a);
% % % [yrot] = sh_rot(yrot, g, b, a);
% % % [zrot] = sh_rot(zrot, g, b, a);
% % % 
% % % %% Step 3: Rotate the ellipsoid representation of the clks too. 
% % % %%%% This is easier ! All we need to do is determine the
% % % %%%% rotation matrix of the ellipsoid in object space (i.e. V above) that takes it to
% % % %%%% canonical form. This is also the rotation matrix that should act on
% % % %%%% the Clks
% % % C = [xrot(:)';yrot(:)';zrot(:)'];
% % % %%%%%%
% % % A = [xrot(2:4); yrot(2:4);zrot(2:4)];
% % % T = A'*A;
% % % [V2 D] = eig(T);
% % % Do = diag(1./diag(sqrt(D)));
% % % cr = Do*A'*C;
% % % X_o = [cr(1,:) cr(2,:) cr(3,:)];

X_o = cs2nocs(X_o);
