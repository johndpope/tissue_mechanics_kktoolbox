function plot_state(u,v,w, F, count, flag, az, el)

cla reset;axis square;warning off;
[t p r] = kk_cart2sph(u,v,w);

%patch('Vertices',[u v w],'Faces',F, 'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');
patch('Vertices',[u v w],'Faces',F, 'FaceVertexCData', t,'FaceColor','flat');

graphlims = [-1.1 1.1];xlim(graphlims);ylim(graphlims); zlim(graphlims);view(3);
if nargin >4, title(num2str(count));view(az,el);end;
hold off;
if nargin >4, zoom(flag);else zoom(2);end
drawnow;warning on;
