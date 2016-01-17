function sh_plot_multiple(xclks, yclks, zclks)
%
%   USAGE:
%   sh_plot_multiple(xclks, yclks, zclks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
gdim = 40;L_max = round(sqrt(length(xclks))-1);
P = partsphere(gdim^2);
x = reshape(P(1,:),gdim,gdim);y = reshape(P(2,:),gdim,gdim);z = reshape(P(3,:),gdim,gdim);[t p r] = kk_cart2sph(x,y,z);
[x y z] = kk_sph2cart(t,p,1);
x = reshape(x,gdim^2,1);y = reshape(y,gdim^2,1);z = reshape(z,gdim^2,1);t = reshape(t,gdim^2,1);p= reshape(p,gdim^2,1);
X = [x y z];
[C] = convhulln(X, {'Qt'});

[L, K ] = indices_gen(1:(L_max + 1)^2); M = length(L);N = length(x);
str = sprintf('gdim_%d_lmax_%d.mat', gdim, L_max);
if exist(str)==2, load(str);disp('Loading precalculated SH basis functions');end;
if ~exist('Y_LK'),
    disp('Calculating basis set functions for the first time.');
    Y_LK  = zeros(N, M);for S = 1:length(L), Y_LK(:,S) = ylk_cos_sin(L(S),K(S),p',t')';end;
    save(str,'Y_LK');
elseif any(size(Y_LK)~=[N M]),
    disp('Recalculating basis set functions to accomodate new truncation/dimension.');
    if size(Y_LK,2)>M && size(Y_LK,1)==N, 
        disp('Truncating precalculated (large) set.');
        Y_LK = Y_LK(:,1:M);
    else,     Y_LK  = zeros(N, M);for S = 1:length(L), Y_LK(:,S) = ylk_cos_sin(L(S),K(S),p',t')';save(str,'Y_LK');end;
    end
else
    disp('No need to recalculate SH basis. Using precalculated set');
end;
% X = [sum(xclks(ones(gdim^2,1),:).*Y_LK,2)  sum(yclks(ones(gdim^2,1),:).*Y_LK,2) sum(zclks(ones(gdim^2,1),:).*Y_LK,2)];
X = Y_LK*[xclks(:) yclks(:) zclks(:)];
%% calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));
Area = sum(twoA)/2; 
F_areas = twoA/2;
n = crossqpr./twoA(:,ones(3,1));
V = sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
disp([Area V v_red]);
Y_LK_ret = Y_LK;
%%%%%% look at shape
figure;patch('Vertices', X, 'Faces', C,'FaceColor', 'red','EdgeColor','none');daspect([1 1 1]);axis off;light; lighting gouraud;lightangle(100,-90)

