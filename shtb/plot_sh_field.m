function plot_sh_field(X_o, sf)
%s = sh_surface(L,sh_basis(L,120));
%s.xc = sh_surface.tr_xc(gtw,L);s = sh_rot(s,0, 0, pi);
gdim = 40;
[xc yc zc] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);

[X, C] = mesh_gen_rand();       % looks better than the icosahedron subdivision
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));

Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
%% generate the shape outline
X = Y_LK(:,1:length(xc))* [xc(:) yc(:) zc(:)];

%% calculate basis for the scalar field
L_max = round(sqrt(length(sf))-1);
[L, K] = indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
Y_LK_sf  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK_sf(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;
sf = Y_LK_sf(:,1:length(sf))*sf(:);%% generate the  scalar field

%%%%%% look at shape
figure;patch('Vertices', X, 'Faces', C,'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
daspect([1 1 1]);
axis off; 
lighting gouraud;
view(69,-43);
%%
function [X, C] = mesh_gen_rand()
dim = 180;
P = partsphere(dim^2);
x = reshape(P(1,:),dim,dim);
y = reshape(P(2,:),dim,dim);
z = reshape(P(3,:),dim,dim);
X = double([x(:) y(:) z(:)]);
[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
