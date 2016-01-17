function [t p] = conformal_spherical_mapping_mc(Xo,F, t, p)
% Xo is the 3-vector of x y z coordinates of the original morphology
% F the connectivity
% t and p the starting values for theta and phi coordinates on sphere

%% configure
gamma = 1e2;      % strain energy modulus --- for area dilation
do = 0.3333;
miu = 1e0;    % shear modulus---- changes in d from do;
newton_step = 1e-5;
niter = 100;
fdiff_step = 0.02;
%% initialize
nvert = length(t);
[x,y,z] = kk_sph2cart(t,p,1);
   %% plot
Xplot = [x(:) y(:) z(:)];
cla;plot_mesh(Xplot,F);view(-60,30);drawnow
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
[Edge, L, face_memb] = edge_info(Xo,F);      % inefficient
disp('Done!');

X = [t(:);p(:)];
Xtry = X;
t(:) = X(1:length(t));
p(:) = X(length(t) + 1:end);
ttry = t;
ptry = p;

for (nix = 1:niter)
    disp(['Iteration: ' num2str(nix)]);
    for cix = 1:length(X)        % loop over the configuration indices (rows of J)
        if mod(cix,1000) == 0,
            disp(['     ' num2str(cix) ' of ' num2str(length(X)) ' finished!']);
        end
        %% the faces that this vertex index is member of.
        if cix<=length(t),  f = face_memb{cix};
        else f = face_memb{cix-length(t)};
        end
        %% calculate the energy before the change
        E = 0;
        for(fix = 1:length(f))  % loop over those faces
            [a d] = tri_spherical_prop([t(F(f(fix),1)) p(F(f(fix),1))], [t(F(f(fix),2)) p(F(f(fix),2))], [t(F(f(fix),3)) p(F(f(fix),3))]);
            E = E + gamma * (a-0.0)^2 + ...
                miu * ((d(1)/sum(d) - do)^2 +(d(2)/sum(d) - do)^2 +(d(3)/sum(d) - do)^2) ;       % energy for entry f(fix) into the column of J (i.e. the triangles face)
        end
        %% apply a change
        Xtry(cix) =  Xtry(cix) + fdiff_step*randn;
        ttry(:) = Xtry(1:length(t)); ptry(:) = Xtry(length(t) + 1:end);
        %% calculate energy after the change
        Etry = 0;
        for(fix = 1:length(f))  % loop over those faces
            [a d] = tri_spherical_prop([ttry(F(f(fix),1)) ptry(F(f(fix),1))], [ttry(F(f(fix),2)) ptry(F(f(fix),2))], [ttry(F(f(fix),3)) ptry(F(f(fix),3))]);
            Etry = Etry + gamma * (a-0.0)^2 + ...
                miu * ((d(1)/sum(d) - do)^2 +(d(2)/sum(d) - do)^2 +(d(3)/sum(d) - do)^2) ;       % energy for entry f(fix) into the column of J (i.e. the triangles face)
            
        end
        if Etry<E, 
            %disp([Etry-E Xtry(cix)-X(cix)]);
            X(cix) = Xtry(cix);%% then make the change
            t(:) = X(1:length(t));
            p(:) = X(length(t) + 1:end);
        end
    end
    
    %% % plot current configuration if desired
    tpl = [ X(1:(nvert))];ppl = [X(nvert+1:end)];ppl = mod(ppl,2*pi);
    u = cos(pi/2-tpl(:)).*cos(ppl(:));v = cos(pi/2-tpl(:)).*sin(ppl(:));w = sin(pi/2-tpl(:));
%     subplot(4,1,1);
    cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', [1 1 1],'FaceColor','flat');axis square;
    colorbar;

    drawnow;

end
t = X(1:length(t));
p = X(length(t) + 1:end);