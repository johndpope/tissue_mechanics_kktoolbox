function [Area,V,v_red,t,p,X,C,Y_LK]=plot_sh2(X1,X2,gdim)
%
%   USAGE:
%   [Area,V,v_red,t,p,X,C,Y_LK]=plot_sh(xclks, yclks, zclks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc = length(X1)/3;yclks = X1(nc+1:2*nc);zclks = X1(2*nc+1:end);xclks = X1(1:nc);
verbose = 1;L_max_x = round(sqrt(length(xclks))-1);L_max_y = round(sqrt(length(yclks))-1);L_max_z = round(sqrt(length(zclks))-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,max([L_max_x L_max_y L_max_z]));
Xx = Y_LK(:,1:length(xclks))*xclks(:);
Xy = Y_LK(:,1:length(yclks))*yclks(:);
Xz = Y_LK(:,1:length(zclks))*zclks(:);
X1 = [Xx Xy Xz];

nc = length(X2)/3;yclks = X2(nc+1:2*nc);zclks = X2(2*nc+1:end);xclks = X2(1:nc);
verbose = 1;L_max_x = round(sqrt(length(xclks))-1);L_max_y = round(sqrt(length(yclks))-1);L_max_z = round(sqrt(length(zclks))-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,max([L_max_x L_max_y L_max_z]));
Xx = Y_LK(:,1:length(xclks))*xclks(:);
Xy = Y_LK(:,1:length(yclks))*yclks(:);
Xz = Y_LK(:,1:length(zclks))*zclks(:);
X2 = [Xx Xy.*1.05 Xz*.9];
%% calculation of the area and volume

%%%%%% look at shape
figure;cla;
patch('Vertices', X1, 'Faces', C,'FaceColor', 'r','FaceAlpha',0.4','EdgeColor','none');
patch('Vertices', X2, 'Faces', C,'FaceColor', 'b','FaceAlpha',1','EdgeColor','none');
daspect([1 1 1]);axis off; lighting gouraud;lightangle(64,-42);lightangle(-100,0);view(69,-43);
















