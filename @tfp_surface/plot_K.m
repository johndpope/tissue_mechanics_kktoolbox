function plot_K(obj)
%% plot on top of the shape the scalar field K (Gaussian curvature)
dfig;
obj = update(obj);
C = reshape(obj.KG,size(obj.x));
surf(double(obj.x),double(obj.y),double(obj.z), C,  'EdgeColor', 'none');
axis equal;axis off;
lighting gouraud;
camlight;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
axis vis3d;


