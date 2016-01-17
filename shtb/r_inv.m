function [X_o X_1 X_2 ang res_p res_o] = r_inv(X_o)
%%% Calculate the Euler angles a,b,g needed for a rotational invariant shape description of X_o shape. 
%%% Input:
%%%     X_o: vector of SHP coefficients in nocs convention.
%%% Output:
%%%     X_o: fully rotationally invariant vector
%%%     X_1: parameterization invariant only (object is not rotated)
%%%     X_2: parameterization and object invariant, but not sign corrected
%%%     ang: transformation that can be used to return X_o into X_1
%%% Uses:
%%%         get_xyz_clks, rotate_x_shp, rotate_z_shp,
%%%         parametric_rotation_objective
%%% Author: Khaled Khairy
%%% To do: Compare this (slow) method with Brechbuehler et al. 1995 method
%%%% WORKS !!! please only modify if you will make a backup copy
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;
mc_iter_p = 500;
plot_flag_p = 0;

mc_iter_o = 500;
plot_flag_o = 0;

[xc yc zc] = get_xyz_clks(X_o);
xpos = xc(1);ypos = yc(1);zpos = zc(1);
[X_1 res_p] = parametrization_invariance(X_o, mc_iter_o, plot_flag_p); if verbose,disp('Canonical parameterization calculated.');end
[X_2 ang res_o] = object_invariance(X_1, mc_iter_o, plot_flag_o);if verbose, disp('Canonical object orientation calculated.');end

[X_o Ihist] = fix_signs(X_2);   % Ihist encodes the inversions needed to reconstruct the original configuration
ang = [ang(:)' Ihist(:)'];
X_o = fix_inversion(X_o);
pass = check_canonicity(X_o);
if pass == 0,
    X_o = fix_canonicity(X_o);
end

%%% put the object at its original position
[xc yc zc] = get_xyz_clks(X_o);
xc(1) = xpos;
yc(1) = ypos;
zc(1) = zpos;

X_o = ([xc(:)' yc(:)' zc(:)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [X_o res_small]= parametrization_invariance(X_o, mc_iter, plot_flag) %#ok<INUSD>
%% parameterization invariance
t = 0;p = 0;nY1 = [ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) ylk_cos_sin(1,0,p, t)/N_LK(1,0) ylk_cos_sin(1,1,p, t)/N_LK(1,1)];
t = pi/2;p = 0;nY2 = [ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) ylk_cos_sin(1,0,p, t)/N_LK(1,0) ylk_cos_sin(1,1,p, t)/N_LK(1,1)];

[xc yc zc] = get_xyz_clks(nocs2cs(tr(X_o,1)));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
res_small = inf;
%%% let's do a quick Monte Carlo estimation to get a good starting value
counter = 0;
while counter<mc_iter
    ang = [rand(1)*2*pi rand(1)*pi rand(1)*2*pi];
    [res] = parametric_rotation_objective(ang,xc,yc,zc, sqrtdD, nY1, nY2,0);
    if res < res_small,res_small = res; ang_min = ang;end
    counter = counter + 1;
end
ang = ang_min;



if res_small > 1e-10
    options =   optimset(...
        'Algorithm', 'levenberg-marquardt','LineSearchType', 'cubicpoly', 'LargeScale', 'off', ...
        'MaxFunEvals', 1e3,'MaxIter', 300,'DiffMaxChange', 0.2,'DiffMinChange', 1e-4,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    
    plot_flag = 0;
    [ang res_small] = fminunc(@parametric_rotation_objective,ang,options,xc,yc,zc, sqrtdD, nY1, nY2,plot_flag);
end
 [xc yc zc] = get_xyz_clks(nocs2cs(X_o));
%%dfig(1);clf;plot_shps((X_o),3,'red');view(3);lighting gouraud;camlight;

[xc] = sh_rot(xc, ang(1), ang(2), ang(3));
[yc] = sh_rot(yc, ang(1), ang(2), ang(3));
[zc] = sh_rot(zc, ang(1), ang(2), ang(3));
X_o = (cs2nocs([xc(:)' yc(:)' zc(:)'])); %  Note the normalization to be compatible with our current basis definition
% dfig(2);clf;plot_shps(X_o,3,'red');view(3);lighting gouraud;camlight;
%%
function [X_out ang res_small] = object_invariance(X_in, mc_iter, plot_flag)
t = 0;p = 0;nY1 = [ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) ylk_cos_sin(1,0,p, t)/N_LK(1,0) ylk_cos_sin(1,1,p, t)/N_LK(1,1)];
t = pi/2;p = 0;nY2 = [ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) ylk_cos_sin(1,0,p, t)/N_LK(1,0) ylk_cos_sin(1,1,p, t)/N_LK(1,1)];

X_in = nocs2cs(X_in);
[xc yc zc] = get_xyz_clks((X_in));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
res_small = inf;
counter = 0;
while counter<mc_iter   %%% let's do a Monte Carlo search to get a  good starting value
    ang = [rand(1)*2*pi rand(1)*pi rand(1)*2*pi];
    res = object_rotation_objective(ang,X_in, sqrtdD,nY1, nY2, plot_flag);
    if res < res_small,res_small = res; ang_min = ang;end
    counter = counter + 1;
end
ang = ang_min;

if res_small > 1e-10    %%% perform nonlinear optimization to refine Monte Carlo guess
    options =   optimset(...
        'Algorithm', 'levenberg-marquardt','LineSearchType', 'cubicpoly', 'LargeScale', 'off', ...
        'MaxFunEvals', 1e3,'MaxIter', 300,'DiffMaxChange', 0.2,'DiffMinChange', 1e-4,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    
    plot_flag = 0;
    [ang res_small] = fminsearch(@object_rotation_objective,ang,options,X_in, sqrtdD, nY1, nY2,plot_flag);
end
ang(1) = mod(ang(1), 2*pi);
ang(2) = mod(ang(2), 2*pi);
ang(3) = mod(ang(3), 2*pi);
X_out = cs2nocs([rotate_shp(X_in,ang)]');
%%
function [X_o Ihist]= fix_signs(X_o)
%%% fix the signs to become canonical (we want the largest values in X to
%%% be -ve, for Y +ve and for Z -ve , as by convention). It also returns
%%% the inversion history so that we can reconstruct the original
%%% orientation from the canonical shape
[xc yc zc] = get_xyz_clks(X_o);
Ihist = [1 1 1];
if sign(xc(4)) == 1, Ihist(1) = -1;xc = -xc;end
if sign(yc(2)) == -1, Ihist(2) = -1;yc = -yc;end
if sign(zc(3)) == 1, Ihist(3) = -1;zc = -zc;end

X_o = [xc(:)' yc(:)' zc(:)'];
%%
function X_o = fix_canonicity(X_o)
%%% if it fails canonicity test, try to fix it
disp('Trying to fix canonicity problem.');
[xc yc zc] = get_xyz_clks(X_o);
if  sign(zc(3))== 1,
    disp('Case 1');
    [xc yc zc] = get_xyz_clks(X_o); X_o = [xc(:)' yc(:)' -zc(:)'];
    pass = check_canonicity(X_o); 
    if pass == 0, 
        [xc yc zc] = get_xyz_clks((X_o)); disp([xc(2:4) yc(2:4) zc(2:4)]);
        warning('Canonical shape transformation inconsistent');
    else
        disp('fixed !');
    end
end

%%
function X_o = fix_inversion(X_o)
%%%% if max(|C20|) is not negative, we invert the shape to make it negative.
%%%% Inversion of SHP shapes occurs by changing the sign of all
%%%% coefficients with even L > 0 coefficients 

 [xc yc zc] = get_xyz_clks(X_o);
 
 vec = ([xc(7) yc(7) zc(7)]);
 indx = find(abs(vec) == max(abs(vec)));    % get the index i (i = x, y, or z) with the largest absolute value
 
 % if this coefficient is -ve, then let us change the sign of all L even
 % coefficients
 if sign(vec(indx))==-1, X_o = fix_xc20(X_o);end




%%













