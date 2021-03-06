function [A, V, Y_LK, all_F,all_X, com, r,dim] = plot_shps(S, nico, option)
%
%   USAGE: S is n x m where n is the number of shapes and m is dimension
%   along which the coefficients are stored (i.e. 3 x that for each
%   coordinate).
%   option: 'red' 'blue' 'volume' 'area' 'random'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, nico = 2; option = 'random';end
if nargin == 2, option = 'rednone';end
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape
L_max = round(sqrt(size(S,2)/3)-1);

% gdim = 20;P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
% [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
% X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

%% using subdivisions of icosahedron
if nico > 6, nico = 6; disp(['Icosahedron subdivision: ' num2str(nico)]);end;

[X,C]=BuildSphere(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
X = double([x(:) y(:) z(:)]);


Y_LK = get_basis(t',p',length(t),max([L_max L_max L_max]));
A =  zeros(size(S,1),1); V = zeros(size(S,1),1);
if nargout>=5
all_F = {};
all_X = {};
com = [];
r = [];
end
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    if isnumeric(option), [a, v,X] = plot_each(xclks, yclks, zclks, Y_LK, C, option(ix));hold on;
    else
        [a, v,X] = plot_each(xclks, yclks, zclks, Y_LK, C, option);hold on;
    end
    A(ix) = a; V(ix) = v;
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
% daspect([1 1 1]);%axis on; 
%axis equal;axis off;lighting flat;%view(2);camlight


dim = length(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area, V,X] = plot_each(xclks, yclks, zclks, Y_LK, C, option)

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
if isnumeric(option),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',option*ones(length(X),1),'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'red'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','k','FaceAlpha',1);
end
if strcmpi(option,'rednone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'rednone_transparent'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',0.4);
end
if strcmpi(option,'greennone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'g','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'greynone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',1);
end

if strcmpi(option,'blue'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'b','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'bluenone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'b','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'volume'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',V,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'area'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',A,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'invisible'),
        % do nothing
end
if strcmpi(option,'random')
     a = rand(1,3);      % random RGB color specification
     if sum(a)> 2.8; a = [1 0 1];end
    patch('Vertices', X, 'Faces', C,'FaceColor', a,'EdgeColor','none','FaceAlpha',1);
end
%lighting phong;


