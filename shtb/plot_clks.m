function plot_clks(clks)
%%% needs to be developed to distinguish different L orders

[xc yc zc] = get_xyz_clks(clks);
figure;clf
subplot(3,1,1); plot(xc, '.');
subplot(3,1,2); plot(yc, '.');
subplot(3,1,3); plot(zc, '.');

disp([xc(:) yc(:) zc(:)]);


% % %%% use icosahedron subdivision
nico = 4;
[X,C]=BuildSphere(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));

% % %%% use partsphere
% % gdim = 15;P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
% % [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
% % X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});



figure;clf;
% % subplot(1,3,1);sh_plot_direct_sh(xc, t,p,X, C);title('x component');view(2);
% % subplot(1,3,2);sh_plot_direct_sh(yc, t,p,X, C);title('y component');view(2);
% % subplot(1,3,3);sh_plot_direct_sh(zc, t,p,X, C);title('z component');view(2);

sh_plot_direct_sh(xc, t,p,X, C);title('x component');view(2);
figure;sh_plot_direct_sh(yc, t,p,X, C);title('y component');view(2);
figure;sh_plot_direct_sh(zc, t,p,X, C);title('z component');view(2);




function sh_plot_direct_sh(c, t,p,X, C)
c = c(:)';
Y_LK = get_basis(t',p',length(t), sqrt(length(c))-1);
r = Y_LK*c';
%r = 1;
% [x y z] = kk_sph2cart(t(:),p(:),r);
% X = [x(:) y(:) z(:)];
patch('Vertices', X, 'Faces', C,'FaceVertexCData',r(:),'FaceColor', 'flat','EdgeColor','k','FaceAlpha',1);
axis equal
view(0,0);
lighting flat;
% camlight
xlabel('x');ylabel('y');zlabel('z');