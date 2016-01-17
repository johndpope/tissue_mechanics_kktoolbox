function [A V v_red] = area_shps(S, nico)
%
%   USAGE: S is n x m where n is the number of shapes and m is dimension
%   along which the coefficients are stored (i.e. 3 x that for each
%   coordinate).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, nico = 4; end
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape
L_max = round(sqrt(size(S,2)/3)-1);
%%%% using subdivisions of icosahedron
if nico > 5, nico = 5;end
[X,C]=BuildSphere(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
Y_LK = get_basis(t',p',length(t),max([L_max L_max L_max]));
A =  zeros(size(S,1),1);
V = zeros(size(S,1),1);
v_red = zeros(size(S,1),1);

for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    [a, v, ~,vr] = calc_each(xclks, yclks, zclks, Y_LK, C);
    A(ix) = a; 
    V(ix) = v;
    v_red(ix) = vr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area, V,X, v_red] = calc_each(xclks, yclks, zclks, Y_LK, C)
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
%% calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);

