function [t p] = conformal_spherical_mapping(Xo,F, t, p)
% Xo is the 3-vector of x y z coordinates of the original morphology
% F the connectivity
% t and p the starting values for theta and phi coordinates on sphere

%% configure
miu = 0.0;    % shear modulus
gamma = 1e2;      % strain energy modulus
newton_step = 1e-5;
niter = 100;
fdiff_step = 0.001;
%% initialize
nvert = length(t);
[x,y,z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)];

%% prepare reference area and edge lengths for every triangle
% we use the triangles of the original morphology for this purpose
ao = zeros(length(F),1);
do = 0.333333 * ones(length(F),3);
%% comment out to use equal and general minimization of area and shear
% % x = Xo(:,1);y = Xo(:,2);z = Xo(:,3);
% % for(fix = 1:length(F))
% %     V1 = [x(F(fix,1)) y(F(fix,1)) z(F(fix,1)) ];
% %     V2 = [x(F(fix,2)) y(F(fix,2)) z(F(fix,2)) ];
% %     V3 = [x(F(fix,3)) y(F(fix,3)) z(F(fix,3)) ];
% %     [a n d] = tri_prop(V1, V2, V3);
% %     [as ds] = tri_spherical_prop([t(F(fix,1)) p(F(fix,1))], [t(F(fix,2)) p(F(fix,2))], [t(F(fix,3)) p(F(fix,3))]);
% %     ao(fix) = a;
% %     do(fix,:)= d/sum(d);        % store the relative triangle length
% % end
% % ao = ao/sum(ao);    % we store the area relative to the total area

                          % multiply by 0.0 to equally minimize all areas
%% calculate preliminary quantities
disp('Gathering edge and link information...');
[Edge, L, face_memb] = edge_info(X,F);      % inefficient
disp('Done!');

%% calculate starting energy
% we use the triangles of the original morphology for this purpose
da = [];
asv = [];
dd = [];
av = [];
E = zeros(length(F),1);
EA = E;
Ed = E;
[x,y,z] = kk_sph2cart(t,p,1);
for(fix = 1:length(F))
    V1 = [x(F(fix,1)) y(F(fix,1)) z(F(fix,1)) ];
    V2 = [x(F(fix,2)) y(F(fix,2)) z(F(fix,2)) ];
    V3 = [x(F(fix,3)) y(F(fix,3)) z(F(fix,3)) ];
%     [a n d] = tri_prop(V1, V2, V3);
    [a d] = tri_spherical_prop([t(F(fix,1)) p(F(fix,1))], [t(F(fix,2)) p(F(fix,2))], [t(F(fix,3)) p(F(fix,3))]);
%     av = [av a];
%     dd = [dd d-ds];
%     da = [da a-as];
%     asv = [asv as];
    
    E(fix) = gamma*(a/4/pi-ao(fix))^2 + ...
        miu * ((d(1)/sum(d) - do((fix),1))^2 +(d(2)/sum(d) - do((fix),2))^2 +(d(3)/sum(d) - do((fix),3))^2) ;       % energy for entry f(fix) into the column of J (i.e. the triangles face)
end
   %% plot
    Xplot = [x(:) y(:) z(:)];
    cla;plot_mesh(Xplot,F);view(-60,30);drawnow
%% start the iterations
Eplus = zeros(length(F),1);
to = t;
po = p;
X = [t(:);p(:)];
J = zeros(length(F),length(X));
teff = [to;to];
peff = [po;po];
for (nix = 1:niter)
    disp(['Iteration: ' num2str(nix)]);
    for cix = 1:length(X)        % loop over the configuration indices (rows of J)
        if mod(cix,1000) == 0,
            disp(['     ' num2str(cix) ' of ' num2str(length(X)) ' finished!']);
        end
        %% apply the change 
%         CHG = sqrt(eps)* X(cix);
        CHG = fdiff_step* X(cix);
        X(cix) = X(cix) + CHG;
        teff(1:length(t)) = X(1:length(t))';
        peff(length(t)+1:end) = mod(X(length(t)+1:end)',2*pi);
        [x,y,z] = kk_sph2cart(teff,peff,1);
         %% calculate the Jacobian values
         % when we change a t or p value there are non-zero entries in J
         % that correspond to the triangels that that configuration index
         % is member of.
        if cix<=length(t),  f = face_memb{cix};
        else f = face_memb{cix-length(t)};
        end
        for(fix = 1:length(f))  % loop over those faces
                V1 = [x(F(f(fix),1)) y(F(f(fix),1)) z(F(f(fix),1)) ];
                V2 = [x(F(f(fix),2)) y(F(f(fix),2)) z(F(f(fix),2)) ];
                V3 = [x(F(f(fix),3)) y(F(f(fix),3)) z(F(f(fix),3)) ];
 %              [a n d] = tri_prop(V1, V2, V3);
                [a d] = tri_spherical_prop([teff(F(f(fix),1)) peff(F(f(fix),1))], [teff(F(f(fix),2)) peff(F(f(fix),2))], [teff(F(f(fix),3)) peff(F(f(fix),3))]);
                Eplus(f(fix)) = gamma * (a-ao(f(fix)))^2 + ...
                                miu * ((d(1)/sum(d) - do(f(fix),1))^2 +(d(2)/sum(d) - do(f(fix),2))^2 +(d(3)/sum(d) - do(f(fix),3))^2) ;       % energy for entry f(fix) into the column of J (i.e. the triangles face)
                EA = gamma * (a-ao(f(fix)))^2;
                %Ed = miu * ((d(1)/sum(d) - do(f(fix),1))^2 +(d(2)/sum(d) - do(f(fix),2))^2 +(d(3)/sum(d) - do(f(fix),3))^2) ;
                J(cix, f(fix)) = (Eplus(f(fix))-E(f(fix)))./CHG;
        end
        X(cix) = X(cix) - CHG;%% restore X to its original value
    end
    %% calculate the solution
    J(isnan(J)) = 0;
    J(isinf(J)) = 100;
    A  = J'*J;b = -E;
    disp(cond(A));
    %% % Calculate the solution
    [dv,~,~,~,~]  = gmres(sparse(A),b,10, [], 200);
    %dv = A\b;
    %[L2,U2] = luinc(sparse(A),1e-6);
    %[dv,flag2,relres2,iter2,resvec2] = bicgstab(A,b,1e-15,10,L2,U2);
    
    disp(norm(b-A*dv)/norm(b));
    
    
    
    %% % and take a step
    step = newton_step.*J*dv;
    %disp(mean(step));
    X = X + newton_step.*J*dv;
    Xcell{nix} = X;
    E(:) = Eplus(:);
    maxE(nix) = max(E(:));
    %% % plot current configuration if desired
    tpl = [ X(1:(nvert))];ppl = [X(nvert+1:end)];ppl = mod(ppl,2*pi);
    u = cos(pi/2-tpl(:)).*cos(ppl(:));v = cos(pi/2-tpl(:)).*sin(ppl(:));w = sin(pi/2-tpl(:));
%     subplot(4,1,1);
    cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', EA,'FaceColor','flat');axis square;
    colorbar;
%     subplot(4,1,2);
%     cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', EA,'FaceColor','flat');axis square;
%     subplot(4,1,3);
%     cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', Ed,'FaceColor','flat');axis square;
%     subplot(4,1,4);
%     cla; plot(maxE);
%     colorbar;
    drawnow;

end
t = X(1:length(t));
p = X(length(t) + 1:end);