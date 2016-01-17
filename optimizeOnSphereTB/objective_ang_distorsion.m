function [R, g] = objective_ang_distorsion(X, flag)
%% Calculates the distorsion of angles from the geodesic square
global F nvert pcount I J gradtheta gradphi temp2 indx beta

if length(X) == 2*nvert,t = X(1:(nvert)); p = X(nvert+1:end); p = mod(p,2*pi);[u, v, w] = kk_sph2cart(t,p ,1);% convert coordinates (on the sphere) to the Cartesian
elseif length(X) == 3*nvert    u = X(1:nvert); v = X(nvert + 1:2*nvert);w = X(2*nvert+1:end);[t,p, r] = kk_cart2sph(u,v,w);end


vert1 = [[u(F(:,1)) v(F(:,1)) w(F(:,1))];[u(F(:,3)) v(F(:,3)) w(F(:,3))]];
vert2 = [[u(F(:,2)) v(F(:,2)) w(F(:,2))];[u(F(:,4)) v(F(:,4)) w(F(:,4))]];
vert3 = [[u(F(:,3)) v(F(:,3)) w(F(:,3))];[u(F(:,1)) v(F(:,1)) w(F(:,1))]];

s1 = vert2-vert1;s2 = vert2-vert3;
ns1 = sqrt(sum(s1.^2,2)); ns2 = sqrt(sum(s2.^2,2));
minvec = abs(sum(s1.*s2./ns1(:,ones(3,1))./ns2(:,ones(3,1)),2));
R = sum(minvec);     % sum the dot products

if nargout>1,   % then calculate the gradient g
       

end
% % Plotting
if flag
    if mod(pcount,flag)==0,
        cla reset;axis square;
        patch('Vertices',[u v w],'Faces',F, 'FaceColor', 'none');
        graphlims = [-1.1 1.1];
        xlim(graphlims);ylim(graphlims); zlim(graphlims);
        % view(3,0);  % good for E
        view(-37,-60);  % good for F
        title(num2str(pcount));
        hold off;
        drawnow;
    end
end
pcount = pcount + 1;

