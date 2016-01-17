function hh = plot_quad_mesh(X,F)

x = X(:,1);
y = X(:,2);
z = X(:,3);
hh = quadmesh(F,x,y,z);
end