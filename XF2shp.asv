function s = XF2shp(X,F, L_max)
if nargin<3, L_max =16;end 
m = surface_mesh(X,F);
m.use_camorbit = 0;
m.newton_step = 0.04;
m.newton_niter = 200;
m.bijective_plot_flag = 1;
%m.optimization_method = 3;
%m = subdivide(m,2);
disp(m);
s = shp_surface(m, L_max); 
