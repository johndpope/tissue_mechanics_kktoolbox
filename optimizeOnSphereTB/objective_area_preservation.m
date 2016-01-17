function [R, g] = objective(X, flag)
%% Calculates the deviation of relative area from that of the original mesh
%% X is of size 2*nvert and corresponds to theta and phi
global F nvert F_areas_o

t = X(1:(nvert)); p = X(nvert+1:end);
[u, v, w] = kk_sph2cart(t,p ,1);% convert coordinates (on the sphere) to the Cartesian

x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2); twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
R = sum(twoA./2/4/pi-F_areas_o);     % sum the dot products

% % Plotting
if flag
    if mod(pcount,flag)==0,
        cla reset;axis square;
        patch('Vertices',[u v w],'Faces',F, 'FaceColor', 'none');
        graphlims = [-1.1 1.1]; xlim(graphlims);ylim(graphlims); zlim(graphlims);
        % view(3,0);  % good for E
        view(-37,-60);  % good for F
        title(num2str(pcount));
        hold off;
        drawnow;
    end
end


