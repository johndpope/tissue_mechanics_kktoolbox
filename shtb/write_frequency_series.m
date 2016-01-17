function write_frequency_series(X_o, gdim, order, filename)
% generates and writes to disk (as tif images) the shapes produced by truncating the X_o series at the
% given orders specified in order
% Usage:
%       write_frequency_series(X_o, gdim, order, filename)
% Example:
%       write_frequency_series(X_o, 120, [1 3 5 22], 'somename')
% Author: Khaled Khairy, January 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4, filename = 'untitled';end
[X,C]=sh_calc(X_o, gdim, order(1));
cla;
fac = 5;
zoomfac = fac-3;
xlim([-max(X(:,1))*fac max(X(:,1))*fac]);
ylim([-max(X(:,2))*fac max(X(:,2))*fac]);
zlim([-max(X(:,3))*fac max(X(:,3))*fac]);

for ix = 1:length(order),
    disp(['Processing Lmax = ' num2str(order(ix))]);
    view([-89 -1]);
    [X,C]=sh_calc(X_o, gdim, order(ix));
    cla;patch('Vertices', X, 'Faces', C,'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    %light('Position',[0.9446 -0.4607 -0.9463]);
    light('Position',[-1.393 0.2456 0]);
    light('Position',[-0.2881 -0.3613 0.8868]);
    daspect([1 1 1]);
    lighting gouraud;
    axis off;
    zoom(zoomfac);
    drawnow;
    order_str = sprintf('%d',order(ix));
    str = [filename '_Lmax_' order_str '.tif'];
    execstr = sprintf('print -dtiff -r600 %s;',str);eval(execstr);
    zoom(1/zoomfac);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,C]=sh_calc(X_o, gdim, Lmax)
X_o = tr(X_o, Lmax);
nc = length(X_o)/3;
xclks = X_o(1:nc);
yclks = X_o(nc+1:2*nc);
zclks = X_o(2*nc+1:end);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,Lmax);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];

