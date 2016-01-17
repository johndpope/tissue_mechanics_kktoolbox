function [t,p,dtline, W] = bijective_map_gen(X, F, L, plotflag, ixN, ixS)
%%% Calculate the bijective mapping of X
if size(X,1) ==3, X = X';end
coord = 1;          % set to 1 for x, 2 for y and 3 for z
if nargin ==4,      % get the indices corresponding to the max and min of z
    z = X(:,coord);
    maxzix = find(z == max(z));minzix = find(z == min(z));
    ixN = maxzix(1);ixS = minzix(1);
%     ixN = 1; ixS = length(X);
end

if plotflag,disp('Calculating theta mapping');end
[t, A, b] = latitude_calc(L, ixN, ixS); %save data_latitude    % Calculate theta (latitude values associated with each vertex)
if plotflag,dfig;patch('Vertices',X,'Faces',F,'FaceVertexCData',t,'FaceColor','interp', 'EdgeColor','k');axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;end


if plotflag,disp('Calculating phi mapping');end
%load data_latitude;
[p, A, b, dtline, W] = longitude_calc(X(:,1), X(:,2), X(:,3),t, A, F, L, ixN, ixS); %save data_longitude

if plotflag
dfig;patch('Vertices',X,'Faces',F,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor','k');axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;
dfig;[u, v, w] = kk_sph2cart(t,p,ones(size(p)));plot_state(u,v,w,F);drawnow
end
% 
