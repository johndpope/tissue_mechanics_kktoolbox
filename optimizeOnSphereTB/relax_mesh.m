function [t,p] = relax_mesh(t, p, E, F, niter, az, el)
%%%%% intitialize the mesh relaxation
nvert = length(t);%disp([length(F) length(E)]);
[u v w] = kk_sph2cart(t,p,1);X = [u v w];
xForces = sparse(nvert,nvert);yForces = sparse(nvert,nvert);zForces = sparse(nvert,nvert);
indx = sub2ind([nvert nvert], E(:,1),E(:,2));
r_length  = sqrt((u(E(:,1))-u(E(:,2))).^2 +(v(E(:,1))-v(E(:,2))).^2 +(w(E(:,1))-w(E(:,2))).^2);
d = zeros(size(r_length));r = [d d d];force_mag = d;Forces = r;TF = zeros(nvert,3);

k = -1;
dt = .01;d_o = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:niter
    d(:) = u(E(:,1)).*u(E(:,2)) + v(E(:,1)).*v(E(:,2)) + w(E(:,1)).*w(E(:,2));d(:) = acos(d);
    r_length(:)  = sqrt((u(E(:,1))-u(E(:,2))).^2 +(v(E(:,1))-v(E(:,2))).^2 +(w(E(:,1))-w(E(:,2))).^2);
    r(:) = [(u(E(:,1))-u(E(:,2))) (v(E(:,1))-v(E(:,2))) (w(E(:,1))-w(E(:,2)))]./r_length(:,ones(3,1));
    force_mag(:) = k*((d-d_o));
    Forces(:) = force_mag(:,ones(3,1)).* r; % the force is calculated here both magnitude(k*(d-d_o)) and direction in the direction of the unit vector r.
    xForces(indx) = Forces(:,1);yForces(indx) = Forces(:,2);zForces(indx) = Forces(:,3);
    TF(:) = [sum(xForces,2) sum(yForces,2) sum(zForces,2)] ;
    X = [u v w];X(2:end,:) = X(2:end,:) + TF(2:end,:).*dt;
    [t, p] = kk_cart2sph(X(:,1), X(:,2), X(:,3));[u v w] = kk_sph2cart(t,p,1);
    if mod(iter,1) == 0;plot_state(u,v,w, F,iter,1,az,el); drawnow;end
end;
