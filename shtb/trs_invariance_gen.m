function [X_o_trans_double_rot_scaled,xrot, yrot,zrot,g,b,a,...
            X_o_not_trans_double_rotated_not_scaled,...
            X_o_not_trans_mono_rotated_not_scaled ] = trs_invariance_gen(clks)
%%%%%% ----------------------------- warning this cannot work !!!
%%% Calculate the Euler angles a,b,g needed for a rotational invariant
%%% shape description of X_o shape.
%%% USAGE:  [X_o,xrot, yrot,zrot,g,b,a ] = rotation_invariance_gen(clks)
xrot = []; yrot = []; zrot = [];X_o = [];
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
c1m1x = xclks(2);  c10x  = xclks(3);  c11x  = xclks(4);
c1m1y = yclks(2);  c10y  = yclks(3);  c11y  = yclks(4);
c1m1z = zclks(2);  c10z  = zclks(3);  c11z  = zclks(4);
Ab = [c1m1x c10x c11x;c1m1y c10y c11y;c1m1z c10z c11z];
T = Ab'*Ab;
[V D] = eig(T);%disp(D);
g = 0;b = 1;a = 0.0; 
ang_o = [g b a]; %starting values
[X] = rotation_objective(ang_o,xclks,yclks,zclks);

options =   optimset(...
            'NonlEqnAlgorithm ', 'gn', 'LineSearchType', 'cubicpoly', ...
            'MaxFunEvals', 5e3,'MaxIter', 30,'DiffMaxChange', .1,'DiffMinChange', 1e-3,...
            'DerivativeCheck','off','LargeScale','off','GradObj','off','GradConstr','off',...
            'Display', 'off','Diagnostics', 'off',...
            'TolCon', 1e-2,'TolFun', 1e-8,'TolX', 1e-10,'MaxSQPIter',100);
% warning off;[ang] = fsolve(@rotation_objective,ang_o,options,xclks,yclks,zclks,[]);warning on;
warning off;[ang] = fsolve(@rotation_objective,ang_o,options,xclks,yclks,zclks);warning on;
g = ang(1);b = ang(2);a = ang(3);%disp(ang);
%  b = pi-b;
%%%%% rotate whole shape
xrot = xclks;yrot = yclks;zrot = zclks;
% xrot = xclks(1:4);yrot = yclks(1:4);zrot = zclks(1:4);
% g1 = -a;b1 = -b;    a1 = -g;
g1 = g;b1 = b;    a1 = a;
[xrot] = sh_rot(xrot, g1, b1, a1);[yrot] = sh_rot(yrot, g1, b1, a1);[zrot] = sh_rot(zrot, g1, b1, a1);
% plot_sh(xrot,yrot,zrot,20);axis on; view(52,23);

%%%% now that we have rotated the clks such that they align with the major
%%%% axes of the ellipsoid, we need to rotate the ellipsoid representation
%%%% of the clks too. This is easier ! All we need to do is determine the
%%%% rotation matrix of the ellipsoid in object space (i.e. V above) that takes it to
%%%% canonical form. This is also the rotation matrix that should act on
%%%% the Clks
% Rg = rot_mx(1,2);
C = [xrot(:)';yrot(:)';zrot(:)'];
%%%%%%
c1m1x = xrot(2);  c10x  = xrot(3);  c11x  = xrot(4);c1m1y = yrot(2);  c10y  = yrot(3);  c11y  = yrot(4);c1m1z = zrot(2);  c10z  = zrot(3);  c11z  = zrot(4);
A = [c1m1x c10x c11x;c1m1y c10y c11y;c1m1z c10z c11z];
T = A'*A;
[V2 D] = eig(T);%disp(D);
Do = diag(1./diag(sqrt(D)));
cr = Do*A'*C;
%%%%%%
scale = max(diag(sqrt(D)));
xrot2 = [0 cr(1,2:end)]./scale;
yrot2 = [0 cr(2,2:end)]./scale;
zrot2 = [0 cr(3,2:end)]./scale;    % translation invariance and scale invariance
fac = 1e2;      % useful when making volume out of this shape
X_o_trans_double_rot_scaled = fac*[xrot2(:)' yrot2(:)' zrot2(:)']';
X_o_not_trans_double_rotated_not_scaled = [cr(1,:) cr(2,:) cr(3,:)];
X_o_not_trans_mono_rotated_not_scaled = [xrot(:)' yrot(:)' zrot(:)'];
% plot_sh(xrot2,yrot2,zrot2,20);axis on; view(52,23);drawnow;
