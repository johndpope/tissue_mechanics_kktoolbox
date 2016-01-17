function [Xvec, F, A, V, Y_LK, C] = export_shp_shapes(S, gdim, filename)
%
%   USAGE:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, gdim = 120; end
verbose = 0;
L_max = round(sqrt(size(S,2)/3)-1);

P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

F = [];
Xvec = [];

Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
A =  []; V = [];
counter = 0;
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    [a, v, X] = calc_each(xclks, yclks, zclks, Y_LK, C, 'r');
    A = [A a]; V = [V v];
    Xvec = [Xvec;X];
    F = [F;C+length(X)*counter];
    counter = counter+1;
end
plot_mesh(Xvec,F);
write_obj([filename '.obj'], Xvec,F);
write_off([filename '.off'], Xvec,F);
write_ply([filename '.ply'], Xvec,F);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area, V, X] = calc_each(xclks, yclks, zclks, Y_LK, C, farbe)
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
%%% calculation of the area and volume
u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);
%%%%%% look at shape
% cla;patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',0.5);
%patch('Vertices', X, 'Faces', C,'FaceVertexCData',V,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
%lighting phong;


