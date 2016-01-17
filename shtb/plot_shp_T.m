function [A, V, Y_LK, C] = plot_shp_T(T, gdim, farbe)
%
%   USAGE:
%   [Area,V,v_red,t,p,X,C,Y_LK]=plot_sh(xclks, yclks, zclks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, gdim = 120; farbe = 'r';end
if nargin == 2, farbe = 'b';end
verbose = 0;
S = T(:,1:(end-1));
L_max = round(sqrt(length(S(1:end)))-1);

P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

 Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
A =  []; V = [];
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
     [a, v] = plot_each(xclks, yclks, zclks, Y_LK, C, T(ix,end));hold on;
%[a, v] = plot_each(xclks, yclks, zclks, Y_LK, C, ix);hold on;
    A = [A a]; V = [V v];
end
hold off
daspect([1 1 1]);axis on; lighting gouraud;
lightangle(64,-42);lightangle(-100,0);view(69,-43);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area, V] = plot_each(xclks, yclks, zclks, Y_LK, C, farbe)

X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];

%% calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);

%%%%%% look at shape
% cla;patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',0.5);
patch('Vertices', X, 'Faces', C,'FaceVertexCData',farbe,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
lighting phong;


