function [Area,V,v_red,t,p,X,C,Y_LK, Xx, Xy, Xz]=plot_sh(xclks, yclks, zclks, gdim, Y_LK)
%
%   USAGE:
%   [Area,V,v_red,t,p,X,C,Y_LK]=plot_sh(xclks, yclks, zclks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,gdim = 40;nc = length(xclks)/3;yclks = xclks(nc+1:2*nc);zclks = xclks(2*nc+1:end);xclks = xclks(1:nc);end
if nargin == 2,gdim = yclks;nc = length(xclks)/3;yclks = xclks(nc+1:2*nc);zclks = xclks(2*nc+1:end);xclks = xclks(1:nc);end
% if nargin ==4,gdim = 40;end
verbose = 0;
L_max_x = round(sqrt(length(xclks))-1);
L_max_y = round(sqrt(length(yclks))-1);
L_max_z = round(sqrt(length(zclks))-1);

P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

if nargin <5, Y_LK = get_basis(t',p',gdim,max([L_max_x L_max_y L_max_z]));end

% clear P x y z t p r 
% Xx = Y_LK(:,1:length(xclks))*xclks(:);
% Xy = Y_LK(:,1:length(yclks))*yclks(:);
% Xz = Y_LK(:,1:length(zclks))*zclks(:);
% X = [Xx Xy Xz];

X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];

%% calculation of the area and volume
% if verbose,
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);
% disp([Area V v_red]);
% end

%%%%%% look at shape
% cla;patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',0.5);
cla;patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','k');
daspect([1 1 1]);axis off; lighting gouraud;lightangle(64,-42);lightangle(-100,0);view(69,-43);
axis on;axis equal

