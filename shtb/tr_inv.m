function [X_o] = tr_inv(X_o)
%%% Calculate the Euler angles a,b,g needed for a L = 1 based rotational invariant
%%% shape description of X_o shape.
[xc yc zc] = get_xyz_clks(X_o);xc(1) = 0;yc(1) = 0;zc(1) = 0;X_o = [xc(:)' yc(:)' zc(:)'];%% translation invariance

% disp('---- Input shape area and volume: ');[A_in V_in] = area_shps_notri(X_o, 30);disp([A_in V_in]);
X_o = parametrization_invariance(X_o);
% disp('---- After parametrization area and volume: ');[A_pout V_pout] = area_shps_notri(X_o, 30);disp([A_pout V_pout]);
% % 
X_o = object_invariance(X_o);
% % % % % 
% % X_o = rearrange_flip_invariance(X_o);
X_o = fix_signs(X_o);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
function X_o = parametrization_invariance(X_o)
%% parametrization invariance
[xc yc zc] = get_xyz_clks(nocs2cs(X_o));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
[X] = parametric_rotation_objective([0 0 0],xc,yc,zc, sqrtdD);

%%% let's do a grid search
% disp('Performing rough grid search ...');
Xsmall = inf;
ang_min = [333 333 333];
g = 0:0.5:2*pi;
b = 0:0.3:pi;
a = 0:0.5:2*pi;
tol = 1e-3;
for aix = 1:length(a)
    %disp([aix/length(a) * 100 ang_min Xsmall]);
    for bix = 1:length(b)
        for gix = 1:length(g)
            ang = [g(gix) b(bix) a(aix)];
            [X] = parametric_rotation_objective(ang,xc,yc,zc, sqrtdD);
            if X < Xsmall,
                Xsmall = X;
                ang_min(:) = ang(:);
                %disp([ang(:)' X]);
            end
            if Xsmall < tol, break;end
        end
        if Xsmall < tol, break;end
    end
    if Xsmall < tol, break;end
end
ang = ang_min;
% disp([ang(:)' Xsmall]);
% disp('Done!   Refining search...');
if Xsmall > 1e-10
    options =   optimset(...
        'NonlEqnAlgorithm ', 'lm', 'LineSearchType', 'cubicpoly', ...
        'MaxFunEvals', 5e3,'MaxIter', 30,'DiffMaxChange', 0.01,'DiffMinChange', 1e-3,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    [ang fval] = fsolve(@parametric_rotation_objective,ang,options,xc,yc,zc, sqrtdD);
    %[ang fval] = fminsearch(@rotation_objective,[0 0 0],options,xc,yc,zc, sqrt([D(2,2)*N_LK(1,0) D(1,1)*N_LK(1,-1)]));
%     disp(['Residual on rotation: ' num2str(fval) ' .']);
end

[xc] = sh_rot(xc, ang(1), ang(2), ang(3));
[yc] = sh_rot(yc, ang(1), ang(2), ang(3));
[zc] = sh_rot(zc, ang(1), ang(2), ang(3));
X_o = (cs2nocs([zc(:)' yc(:)' xc(:)'])); %  Note the normalization to be compatible with our current basis definition
%%
function X_out = object_invariance(X_in)
[xc yc zc] = get_xyz_clks(nocs2cs(X_in));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
[X] = object_rotation_objective([0 0 0],X_in, sqrtdD);
%%% let's do a grid search
% disp('Performing rough grid search ...');
Xsmall = inf;
ang_min = [333 333 333];
g = 0:0.5:2*pi;
b = 0:0.3:pi;
a = 0:0.5:2*pi;
tol = 1e-3;
for aix = 1:length(a)
    %disp([aix/length(a) * 100 ang_min Xsmall]);
    for bix = 1:length(b)
        for gix = 1:length(g)
            ang = [g(gix) b(bix) a(aix)];
            [X] = object_rotation_objective(ang,X_in, sqrtdD);
            if X < Xsmall,
                Xsmall = X;
                ang_min(:) = ang(:);
                %disp([ang(:)' X]);
            end
            if Xsmall < tol, break;end
        end
        if Xsmall < tol, break;end
    end
    if Xsmall < tol, break;end
end
ang = ang_min;
% disp([ang(:)' Xsmall]);
% disp('Object space rotation:   Done!   Refining search...');
if Xsmall > 1e-10
    options =   optimset(...
        'NonlEqnAlgorithm ', 'lm', 'LineSearchType', 'cubicpoly', ...
        'MaxFunEvals', 5e3,'MaxIter', 30,'DiffMaxChange', 0.01,'DiffMinChange', 1e-3,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    [ang fval] = fsolve(@object_rotation_objective,ang,options,X_in, sqrtdD);
    %[ang fval] = fminsearch(@rotation_objective,[0 0 0],options,xc,yc,zc, sqrt([D(2,2)*N_LK(1,0) D(1,1)*N_LK(1,-1)]));
% disp(['Object space rotation:   Residual on rotation: ' num2str(fval) ' .']);
end

% [X] = object_rotation_objective(ang,X_in, sqrtdD);

X_out = rotate_shp(X_in,ang);

% % [xc yc zc] = get_xyz_clks(times_Basis(X_in));    % negate the normalization
% % C = [xc(:)'; yc(:)'; zc(:)'];
% % Ar = [xc(2:4)';yc(2:4)';zc(2:4)'];
% % [V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
% % Cr = V'*C;  % by applying the transformation to allcoefficients
% % 
% % 
% % % dd = sqrt(diag(D));
% % % D = diag(D);D = diag(sqrt(D));Rxyz = D.*Ar'*Ar;
% % % Cr = V'*[1 1 1]';
% % %
% % xc = Cr(1,:);yc = Cr(2,:);zc = Cr(3,:);
% % X_out = over_Basis([xc(:)' yc(:)' zc(:)']); %(weighting by the basis normalization)
% % plot_sh_notri(X_out);
%%
function X_o = rearrange_flip_invariance(X_o)
%% %% Rearrange and flip signs
[r_a c_a] = check_position(X_o,1);
[r_b c_b] = check_position(X_o,2);
[r_g c_g] = check_position(X_o,3);
%disp([r_a c_a;r_b c_b;r_g c_g]);
% % % we need to permute and rotate to obtain the largest value as the upper left
% % % ellipsoid parametric equation matrix element.
% % % put the largest value in the upper left
% place the largest onto the first row
if r_a == 2
    [xc yc zc] = get_xyz_clks(X_o);
    X_o = [yc(:)' xc(:)' zc(:)'];
elseif r_a==3
    [xc yc zc] = get_xyz_clks(X_o);
    X_o = [zc(:)' xc(:)' yc(:)'];
end
% we now have to rotate such that we get largest element onto the top
% left element, i.e. element number 1.
if     c_a ==2,
    X_o = rotate_parameters(X_o, 3, pi/2);
    X_o = fix_2by2(X_o);
elseif c_a ==3,
    X_o = rotate_parameters(X_o, 3, pi/2);      %%% this one still needs work
    X_o = fix_2by2(X_o);
end
%%
function X_o = fix_signs(X_o)
%%% fix the signs 
[xc yc zc] = get_xyz_clks(X_o);
xc = sign(xc(find(abs(xc)==max(abs(xc)))))*xc(:);
yc = sign(yc(find(abs(yc)==max(abs(yc)))))*yc(:);
zc = sign(zc(find(abs(zc)==max(abs(zc)))))*zc(:);
X_o = [xc(:)' yc(:)' zc(:)'];
%%
function X_o = rotate_parameters(X_in, ax , angle)

if ax == 1,      % rotate around x
    [xc yc zc] = get_xyz_clks(nocs2cs(X_in));
    a = 0;    % rotate around y(third rotation)
    b = 0;    % rotate around x(second rotation)
    g = angle;    % rotate around y(first rotation)
    xc = sh_rot(xc,g,b,a);yc = sh_rot(yc,g,b,a);zc = sh_rot(zc,g,b,a);
    X_o = cs2nocs([xc(:)' yc(:)' zc(:)']);

elseif ax == 2,      % rotate around y
    [xc yc zc] = get_xyz_clks(nocs2cs(X_in));
    a = 0;    % rotate around y(third rotation)
    b = angle;    % rotate around x(second rotation)
    g = 0;    % rotate around y(first rotation)
    xc = sh_rot(xc,g,b,a);yc = sh_rot(yc,g,b,a);zc = sh_rot(zc,g,b,a);
    X_o = cs2nocs([xc(:)' yc(:)' zc(:)']);

elseif ax == 3,      % rotate around z
    [xc yc zc] = get_xyz_clks(nocs2cs(X_in));
    a = pi/2;    % rotate around y(third rotation)
    b = angle;    % rotate around x(second rotation)
    g = pi/2;    % rotate around y(first rotation)
    xc = sh_rot(xc,g,b,a);yc = sh_rot(yc,g,b,a);zc = sh_rot(zc,g,b,a);
    X_o = cs2nocs([xc(:)' yc(:)' zc(:)']);
end
%%
function [r c] = check_position(X_o,pos)
[xc yc zc] = get_xyz_clks(times_Basis(X_o));

xc = sign(xc(find(abs(xc)==max(abs(xc)))))*xc(:);
yc = sign(yc(find(abs(yc)==max(abs(yc)))))*yc(:);
zc = sign(zc(find(abs(zc)==max(abs(zc)))))*zc(:);
Arp = [xc(2:4)'; yc(2:4)'; zc(2:4)'];

if pos == 1,        % check the alpha position
[r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the largest value ?
elseif pos ==2,     % check the beta position
    [r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the largest value ?
    Arp(:,c) = 0;
    [r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the second largest value ?
elseif pos ==3,     % check the beta position
    [r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the largest value ?
    Arp(:,c) = 0;
    [r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the second largest value ?
    Arp(:,c) = 0;
    [r c] = ind2sub(size(Arp), find(Arp(:)==max(Arp(:)))); %which column has the third largest value ?
end
%%
function X_o = fix_2by2(X_o)
[r_b c_b] = check_position(X_o,2);
if r_b==3,
    [xc yc zc] = get_xyz_clks(X_o);
    X_o = [xc(:)' zc(:)' yc(:)'];
end

if c_b ==3,
    X_o = rotate_parameters(X_o, 3, pi/2);
end
%%



