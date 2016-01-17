function [Area,V,v_red,t,p,X,C,Y_LK]=shape_explorer_plot_sh(X_o, Y_LK, C, axis1)
%
%   USAGE:
%   [Area,V,v_red,t,p,X,C,Y_LK]=plot_sh(xclks, yclks, zclks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xclks yclks zclks] = shp_surface.get_xyz_clks(X_o);
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
%% calculation of the area and volume

u = X(:,1); v = X(:,2); w = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2; 
F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
V = abs(V);v_red = abs(v_red);
% disp([Area V v_red]);


%%%%%% look at shape
cla;axes(axis1);patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','k');daspect([1 1 1]);axis off;light; lighting gouraud;lightangle(100,-90)
rotate3d on;
