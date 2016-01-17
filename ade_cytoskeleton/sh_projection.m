function [A, V, h, v, xclks, yclks, zclks, C, F_areas] = sh_projection(gdim, xL_max, yL_max, zL_max, x, y, z, t, p, plotflag, shape_prop)
%%% let us calculate the shape properties using the multiple expansion
%%% This involves the expansion of the three functions x(t,p), y(t,p) and z(t,p) on the sphere.
%%% and calculation of properties based on this expansion.
%%% We will assume that the vector field to be expanded is contained in
%%% [xyz] and the vectors are defined on the sphere at the spherical
%%% coordinated t and p.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 0; V = 0; h = 0; v = 0;C = []; F_areas = [];
x = x - mean(x);y = y - mean(y);z = z -mean(z);

disp('Spherical harmonics projection');
A = 0; V = 0; h = 0; v = 0;verbose = 0;
%%%  Expansion into X(t,p), Y(t,p) and Z(t,p)
L_max   = max([xL_max yL_max zL_max]);
% [xclks] = sh_expand(x, xL_max,t,p);if verbose,disp('x expanded');end
% [yclks] = sh_expand(y, yL_max,t,p);if verbose,disp('y expanded');end
% [zclks] = sh_expand(z, zL_max,t,p);if verbose,disp('z expanded');end
[L, K ] = indices_gen(1:(L_max + 1)^2); 
M = length(L); %% number of functions in expansion
N = length(x); %% number of data points
A  = zeros(N, M);
for S = 1:length(L), A(:,S) = ylk_cos_sin_old(L(S),K(S),p,t)';end;warning off MATLAB:divideByZero;
[U, S, V] = svd(A, 0);invS = 1./(S);invS(invS==inf) = 0;
xclks = (V*invS) * (U'*x');yclks = (V*invS) * (U'*y');zclks = (V*invS) * (U'*z');

xclks((L_max + 1)^2 + 1)  = 0; xclks = xclks(1:end-1);
yclks((L_max + 1)^2 + 1)  = 0; yclks = yclks(1:end-1);
zclks((L_max + 1)^2 + 1)  = 0; zclks = zclks(1:end-1);
if shape_prop || plotflag
%%% preparation of quantities needed for the calculation of the shape properties
global Y_LK P_LK Y_LK_phi Y_LK_theta P_LK_T Y_PP Y_TT Y_TP wp wt
[t wt]                  = gaussquad(gdim, 0, pi);
[p wp]                  = gaussquad(gdim,0,2*pi);
[p t]                   = meshgrid(p,t);
[wp wt]                 = meshgrid(wp, wt);
[Y_LK]                  = precalc_ylk_cos_sin_old(p, t, L_max);


c = xclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
x = sum(c.*Y_LK,3);%xp =(sum(c.*Y_LK_phi,3));xt = (sum(c.*Y_LK_theta,3));xpp = (sum(c.*Y_PP,3));xtt = (sum(c.*Y_TT,3));xtp = (sum(c.*Y_TP,3));
c = yclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
y = sum(c.*Y_LK,3);%yp =(sum(c.*Y_LK_phi,3));yt = (sum(c.*Y_LK_theta,3));ypp = (sum(c.*Y_PP,3));ytt = (sum(c.*Y_TT,3));ytp = (sum(c.*Y_TP,3));
c = zclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
z = sum(c.*Y_LK,3);%zp =(sum(c.*Y_LK_phi,3));zt = (sum(c.*Y_LK_theta,3));zpp = (sum(c.*Y_PP,3));ztt = (sum(c.*Y_TT,3));ztp = (sum(c.*Y_TP,3));

if shape_prop, 
    X = [reshape(x,gdim^2,1), reshape(y,gdim^2,1),reshape(z,gdim^2,1)];         % three-vector of resampled coordinates
    
    [u, v, w] = kk_sph2cart(t,p , 1);% convert coordinates (on the sphere) to the Cartesian
    Y = [reshape(u, length(X), 1) reshape(v, length(X), 1) reshape(w, length(X), 1)];
    C = convhulln(Y);%% Triangulate
    H = [];
    [A, V, v, F_areas, h, H, Eb, da] = triangulated_props(X, C, 1);
%     [A, V, v, F_areas] = triangulated_props(X, C, 1);
    if ~isempty(H), H = reshape(H,gdim,gdim);end;
     disp([A V v h da Eb]);
else C = [];
end

%%%%%% look at shape
if plotflag
    if     plotflag ==1 surf(x,y,z,sqrt(x.^2+y.^2+z.^2));grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    elseif plotflag ==2, surf(x,y,z,t);grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    elseif plotflag ==3, surf(x,y,z,p);grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    end
end
end
