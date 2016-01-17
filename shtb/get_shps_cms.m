function [all_F,all_X, com, r] = get_shps_cms(S)
%
%   USAGE: S is n x m where n is the number of shapes and m is dimension
%   along which the coefficients are stored (i.e. 3 x that for each
%   coordinate).
%   option: 'red' 'blue' 'volume' 'area' 'random'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nico = 3;
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape
L_max = round(sqrt(size(S,2)/3)-1);
%%%% using subdivisions of icosahedron
[X,C]=BuildSphere(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
Y_LK = get_basis(t',p',length(t),max([L_max L_max L_max]));
all_F = {};
all_X = {};
com = [];
r = [];
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
    all_F{ix} = C;
    all_X{ix} = X;
    c = sum(X)./length(X);
    com = [com;c];
    distances = sqrt(sum( (X-c(ones(size(X,1),1),:)).^2,2));
    r  = [r;mean(distances)];
end

