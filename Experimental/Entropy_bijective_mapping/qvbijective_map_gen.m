function [t,p, E, L, face_memb, A,b,Aq, bq, dl] = qvbijective_map_gen(X,F, ixN, ixS)
%%% Calculate the bijective mapping
[E,L,face_memb, Ld, Li] = edge_info(X,F);
[t, A, b] = qvlatitude_calc(Ld, ixN, ixS); %save data_latitude    % Calculate theta (latitude values associated with each vertex)
plot_result(X,F,t);
[p, Aq, bq, dl] = qvlongitude_calc_02(X,t, A, F, L, ixN, ixS, face_memb, Ld, Li); %save data_longitude
figure;plot_result(X,F,p);
figure;plot_result_on_sphere(t,p,F);
%%
function plot_result(X,F,a)
% h = dfig;
clf;
patch(  'vertices',X(:,1:3),'faces',F,...
            'FaceColor','interp', 'FaceVertexCData',a,...
            'CDataMapping','scaled',...
            'EdgeColor','k');
axis equal;axis off;
% set(h,'Color','Black');
% camlight;view(90,-90);camlight;colorbar
% lighting flat;cameramenu
%%
function plot_result_on_sphere(t,p,F)
% h = dfig;
[u v w] = kk_sph2cart(t,p,1);

clf;
patch(  'vertices',[u v w],'faces',F,...
            'FaceColor','interp', 'FaceVertexCData',t,...
            'CDataMapping','scaled',...
            'EdgeColor','k');
axis equal;axis off;
% set(h,'Color','Black');
% camlight;view(90,-90);camlight;colorbar
% lighting flat;cameramenu