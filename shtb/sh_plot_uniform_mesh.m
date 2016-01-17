function [xr, yr, zr, C] = sh_plot_uniform_mesh(c, gdim)
%%% just plot the radial function given its SH coefficients
%%% i.e. we are performing spherical harmonic synthesis
if nargin == 1, gdim = 80;end
L_max = (sqrt(length(c))-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,L_max);

r = Y_LK(:,1:length(c))*c(:);
xr = r.*sin(t(:)).*cos(p(:));
yr = r.*sin(t(:)).*sin(p(:));
zr = r.*cos(t(:));
% figure(1);cla;
patch('Vertices', [xr(:) yr(:) zr(:)],'FaceVertexCData', zr(:),...
      'Faces', C,'FaceColor', 'interp','EdgeColor','k',...
      'FaceAlpha',.3);
daspect([1 1 1]);axis off;light; lighting gouraud;lightangle(100,-90);view(3);
