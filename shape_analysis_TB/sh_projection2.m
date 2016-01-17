function [Area, V, h, v, xclks, yclks, zclks, X, C,t,p, A] = ...
    sh_projection2(gdim, xL_max, yL_max, zL_max, x, y, z, t, p)
%%% The expansion of the three functions x(t,p), y(t,p) and z(t,p) on the sphere.
%%% [xyz] and the vectors are defined on the sphere at the spherical
%%% coordinates t and p.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent Y_LK
%if nargin<10, verbose = 1;end
verbose = 0;
%x = x - mean(x);y = y - mean(y);z = z -mean(z);

if verbose, disp('Spherical harmonics projection');end
Area = 0; V = 0; h = 0; v = 0;
L_max   = xL_max;
[L, K ] = indices_gen(1:(L_max + 1)^2); 
M = length(L); %% number of functions in expansion
N = length(x); %% number of data points
% disp('Setting up basis functions matrix');
A  = zeros(N, M, 'single');
if verbose, h = waitbar(0,'Setting up basis functions matrix...');end
for S = 1:length(L),
%     if ~mod(S,50),disp([num2str(S) ' of ' num2str(length(L))]);end;
%     A(:,S) = ylk_cos_sin_nocs(L(S),K(S),p(:)',t(:)')';  % use nocs version
    A(:,S) = ylk_cos_sin(L(S),K(S),p(:)',t(:)')';  % use nocs version
    if verbose, waitbar(S/length(L),h);end
end;
if verbose, close(h);end
if verbose, save data_temp_A A x y z L_max gdim;end
h = 0;
warning off MATLAB:divideByZero;
if verbose, disp('Solving linear system by singular value decomposition');end
[U, S, V] = svd(A, 'econ');invS = 1./(S);invS(invS==inf) = 0;
xclks = (V*invS) * (U'*x(:));yclks = (V*invS) * (U'*y(:));zclks = (V*invS) * (U'*z(:));
if verbose, save data_temp_A  x y z gdim xclks yclks zclks L_max;end

if verbose, disp('Plotting shape based on fitted parameters (CLKs)');end
X_o = [xclks(:)' yclks(:)' zclks(:)'];
% X_o = cs2nocs(X_o);
nc = round(length(X_o)/3);xclks = X_o(1:nc);yclks = X_o(nc+1:2*nc); zclks = X_o(2*nc+1:3*nc);
%plot_sh_notri(X_o);

% figure;X_o = [xclks(:);yclks(:);zclks(:)];
if verbose, [Area,V,v,t,p,X,C,Y_LK_ret]=plot_sh(X_o);
else
    [xclks yclks zclks] = get_xyz_clks(X_o);
    L_max_x = round(sqrt(length(xclks))-1);
    L_max_y = round(sqrt(length(yclks))-1);
    L_max_z = round(sqrt(length(zclks))-1);
    P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
    [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
    X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
    Y_LK = get_basis(t',p',gdim,max([L_max_x L_max_y L_max_z]));
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
    %% calculation of the area and volume
    u = X(:,1); v = X(:,2); w = X(:,3);
    crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
    twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2;
    F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
    V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
    Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
    V = abs(V);v_red = abs(v_red);

end

