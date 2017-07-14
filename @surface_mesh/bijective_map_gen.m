function [t,p,dtline, W, A, b] = bijective_map_gen(X, F, L, plotflag, ixN, ixS)
%%% Calculate the bijective mapping of X
if size(X,1) ==3, X = X';end
if nargin <5      % get the indices corresponding to the max and min of z
    for coord = 1:3          % set to 1 for x, 2 for y and 3 for z
        c = X(:,coord);
        maxcix = find(c == max(c));
        mincix = find(c == min(c));
        cixN(coord) = maxcix(1);    % store the proposed Northpole
        cixS(coord) = mincix(1);    % store the proposed Southpole
        d(coord) = abs(c(cixN(coord))-c(cixS(coord)));    % store the distance between the two
    end
    %indx = find(d == max(d));indx = indx(1);
indx = 3;
    ixN = cixN(indx);
    ixS = cixS(indx);
end

%if plotflag,disp('Calculating theta mapping');end
[t, A, b] = surface_mesh.latitude_calc(L, ixN, ixS); %save data_latitude    % Calculate theta (latitude values associated with each vertex)
if plotflag, figure;patch('Vertices',X,'Faces',F,'FaceVertexCData',t,'FaceColor','interp', 'EdgeColor','k');axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;end


%if plotflag,disp('Calculating phi mapping');end
%load data_latitude;
[p, A, b, dtline, W] = surface_mesh.longitude_calc(X(:,1), X(:,2), X(:,3),t, A, F, L, ixN, ixS); %save data_longitude

if plotflag
figure;patch('Vertices',X,'Faces',F,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor','k');axis square;daspect([1 1 1]);rotate3d;view(3);drawnow;
figure;[u, v, w] = kk_sph2cart(t,p,ones(size(p)));
surface_mesh.plot_state(u,v,w,F);
figure;hist(t,100);title('theta histogram');
figure;hist(p,100);title('phi histogram');
drawnow

end
% 
