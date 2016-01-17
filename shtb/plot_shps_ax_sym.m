function [A, V, Y_LK, all_F,all_X, com, r] = plot_shps_ax_sym(S, gdim, option)
%
%   USAGE: S is n x m where n is the number of shapes and m is dimension
%   along which the coefficients are stored (i.e. 3 x that for each
%   coordinate).
%   option: 'red' 'blue' 'volume' 'area' 'random'
% Input: S is assumed to be related to the normal xclks, yclks and zclks by
%        S = [xclks(4) yclks(2) zclks(k0_ixs)]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, gdim = 40; option = 'red';end
if nargin == 2, option = 'red';end
verbose = 0;
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape

L_max = length(S(3:end));

P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));

%%% choos the basis with axial symmetry
L_list = [1:L_max];                 % let's neglect L = 0;
k0_ixs = L_list.^2 + L_list + 1;    % indices with k = 0;
Y_LK_xy = Y_LK(:,1:4);              % for xy we only need up to L = 1 coefficients specifically for x we need L = 1K = 1, for y we need L = 1 K = -1
Y_LK = Y_LK(:,k0_ixs);


A =  []; V = [];
if nargout>=5
all_F = {};
all_X = {};
com = [];
r = [];
end
for ix = 1:size(S,1), %loop over the shapes
    %[xclks yclks zclks] = get_xyz_clks(S(ix,:));
    xclks = S(1);   % L = 1 K = 1;
    yclks = S(2);   % L = 1 K = -1
    zclks = S(3:end);   % from L = 1, K = 0 to all other K = 0 ones
    [a, v,X] = plot_each(xclks, yclks, zclks, Y_LK_xy, Y_LK, C, option);hold on;
    A = [A a]; V = [V v];
    if nargout>=5,
       all_F{ix} = C;
       all_X{ix} = X;
       c = sum(X)./length(X);
       com = [com;c];
       
       distances = sqrt(sum( (X-c(ones(size(X,1),1),:)).^2,2));
       r  = [r;mean(distances)];
    end
end
hold off
daspect([1 1 1]);axis on; 
%lighting gouraud;
%lightangle(64,-42);lightangle(-100,0);view(69,-43);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area, V,X] = plot_each(xclks, yclks, zclks, Y_LK_xy, Y_LK, C, option)

x = Y_LK_xy(:,4)*xclks;
y = Y_LK_xy(:,2)*yclks;
z = Y_LK*zclks(:);
X = [x(:) y(:) z(:)];

%% calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);
cla
%%%%%% look at shape
if isnumeric(option),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',option,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'red')||strcmpi(option,'r'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','k','FaceAlpha',1);
end
if strcmpi(option,'blue')||strcmpi(option,'b'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'b','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'volume'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',V,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'area'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',A,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'random')
     a = rand(1,3);      % random RGB color specification
     if sum(a)> 2.8; a = [1 0 1];end
    patch('Vertices', X, 'Faces', C,'FaceColor', a,'EdgeColor','none','FaceAlpha',1);
end
%lighting phong;


