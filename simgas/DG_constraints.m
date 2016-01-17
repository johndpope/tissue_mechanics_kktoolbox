function [c, ceq] = DG_constraints(X_o, phys)

verbose = 0;
A_o = phys.A_o;
V_o = phys.V_o;
% a_o_bar = phys.a_o_bar;
v_red_o = phys.v_red_o;

Y = phys.B.Y;
Y_P = phys.B.Y_P;
Y_T = phys.B.Y_T;
Y_PP = phys.B.Y_PP;
Y_TT = phys.B.Y_TT;
Y_TP = phys.B.Y_TP;
p    = phys.B.p;
t    = phys.B.t;
w    = phys.B.w;

%% 
X_o = recover_parms(X_o,phys);
X_o = over_Basis(X_o);
[xclks yclks zclks] = get_xyz_clks(X_o);
X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
[res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);% self-intersection test
%res_self = 0;
flag = 1;
[res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm, flag);% master intersection test
% res_self=0; res_vm = 0;
%%
if res_self == 0 && res_vm == 0
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gdimp = length(p);gdimt = length(t);
    c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    x = sum(c.*Y,3);
    xp =(sum(c.*Y_P,3));
    xt = (sum(c.*Y_T,3));
    xpp = (sum(c.*Y_PP,3));
    xtt = (sum(c.*Y_TT,3));
    xtp = (sum(c.*Y_TP,3));
    
    c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    y = sum(c.*Y,3);
    yp =(sum(c.*Y_P,3));
    yt = (sum(c.*Y_T,3));
    ypp = (sum(c.*Y_PP,3));
    ytt = (sum(c.*Y_TT,3));
    ytp = (sum(c.*Y_TP,3));
    
    c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    z = sum(c.*Y,3);
    zp =(sum(c.*Y_P,3));
    zt = (sum(c.*Y_T,3));
    zpp = (sum(c.*Y_PP,3));
    ztt = (sum(c.*Y_TT,3));
    ztp = (sum(c.*Y_TP,3));
    %%%%%%%%%%%%%%%%%%%%%%%% Calculate the geometrically important quantities
    
    X  =[x(:) y(:) z(:)];
    Xt = [xt(:) yt(:) zt(:)];
    Xp = [xp(:) yp(:) zp(:)];
    Xpp = [xpp(:) ypp(:) zpp(:)];
    Xtp = [xtp(:) ytp(:) ztp(:)];
    Xtt = [xtt(:) ytt(:) ztt(:)];
    
    E = dot(Xt,Xt,2);
    F = dot(Xt,Xp,2);
    G = dot(Xp,Xp,2);
    SS = (cross(Xt,Xp,2));
    SSn = sqrt(E.*G-F.*F);
    n = SS./SSn(:,ones(1,3));
%     L = dot(Xtt,n,2);
%     M = dot(Xtp,n,2);
%     N = dot(Xpp,n,2);
    
    %% calculation of geometric properties
    V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
    A = sum(w.*SSn);
    % % % H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
    % % % K = (L.*N - M.*M)./(E.*G-F.*F);
    % % % T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)
    % % %
    % % % M = sum(sum(w.*H(:).*SSn));         % as in Seifert et al. '91
    % % % a_bar = -M/(4*pi*sqrt(A/4/pi));
    
    Vo = 4/3*pi*(A/4/pi)^(3/2);
%     v_red = V/Vo;
    %%
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [];
    ix = 1;
     ceq(ix) = (1-V/V_o); ix = ix + 1; % volume constraint
%     ceq(ix) = (1-v_red/v_red_o); ix = ix + 1;
%     ceq(ix) = (1-A/A_o); ix = ix + 1;
    ceq(ix) = 0;ix = ix + 1;
else
    if verbose, disp(['Constraint violation: ' num2str([res_self res_vm any(zviolation)])]);end
    %% quick calculation of the area and volume using the triangular mesh
    u = X(:,1); v = X(:,2); w = X(:,3);
    C = phys.C_slave;
    x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
    q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
    crossqpr = cross(q,r,2);
    twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
    A = sum(twoA)/2;                    % this is the total area
%     F_areas = twoA/2;                   % this is the vector of face areas
    n = crossqpr./twoA(:,ones(3,1));
    V = -sum(1/3*(dot(n,[x1 y1 z1], 2).*twoA./2));
    Vo = 4/3*pi*(A/4/pi)^(3/2);
%     v_red = V/Vo;
    %%
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [];
    ix = 1;
     ceq(ix) = (1-V/V_o); ix = ix + 1;
%     ceq(ix) = (1-v_red/v_red_o); ix = ix + 1;
    ceq(ix) = res_vm *1e6;ix = ix + 1;
end
