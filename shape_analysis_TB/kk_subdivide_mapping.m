function [X1,tn,pn, faces] = kk_subdivide_mapping(X,t, p, F, ixlist)
%% subdivide the mapping on the sphere
if nargin<5, ixlist = 1:size(F,1);end
p = mod(p,2*pi);
faces = [];
if size(X,2)== 3
    X1 = X;
else
    X1 = X';
end
tn = t(:);
pn = p(:);
    % we need to subdivide only the given triangles
    faces = F;
    for ix = ixlist(:)',
        newV = [(X(F(ix,1),1)+  X(F(ix,2),1) + X(F(ix,3),1))/3 ...
            (X(F(ix,1),2) + X(F(ix,2),2) + X(F(ix,3),2))/3 ...
            (X(F(ix,1),3) + X(F(ix,2),3) + X(F(ix,3),3))/3];
        X1 = [X1;newV];
        
        %%% determine triangle on the sphere
        [u v w] = kk_sph2cart(...
            [t(F(ix,1)); t(F(ix,2)); t(F(ix,3))],...
            [p(F(ix,1)); p(F(ix,2)); p(F(ix,3))],...
            1);
        SV = [sum(u) sum(v) sum(w)]/3;  % calculate center of mass
        [newt newp r] = kk_cart2sph(SV(1), SV(2), SV(3));
        tn = [tn;newt];
        pn = [pn;newp];
        newix = size(X1,1);
        f1 = [F(ix,1) F(ix,2) newix];
        f2 = [F(ix,2) F(ix,3) newix];
        f3 = [F(ix,3) F(ix,1) newix];
        faces = [faces;f1;f2;f3];
        
% %         u = cos(pi/2-tn(:)).*cos(pn(:));v = cos(pi/2-tn(:)).*sin(pn(:));w = sin(pi/2-tn(:));
% %         U = [u(:) v(:) w(:)];
% %         clf; patch('Vertices',U,'Faces',[f1;f2;f3],'FaceVertexCData', t,...
% %     'FaceColor','flat');axis square;

    end
    faces(ixlist,:) = [];
    

pn = mod(pn,2*pi);
% % figure;
% % u = cos(pi/2-t(:)).*cos(p(:));v = cos(pi/2-t(:)).*sin(p(:));w = sin(pi/2-t(:));
% % U = [u(:) v(:) w(:)];
% % clf; patch('Vertices',U,'Faces',F,'FaceVertexCData', t,...
% %     'FaceColor','flat');axis square;
% % view(-70,4);colorbar;
% % zoom(1);drawnow;title('old spherical topology');
% % 
% % figure;
% % u = cos(pi/2-tn(:)).*cos(pn(:));v = cos(pi/2-tn(:)).*sin(pn(:));w = sin(pi/2-tn(:));
% % U = [u(:) v(:) w(:)];
% % clf; patch('Vertices',U,'Faces',faces,'FaceVertexCData', tn,...
% %     'FaceColor','flat');axis square;
% % view(-70,4);colorbar;
% % zoom(1);drawnow;title('New spherical topology');
% % figure;
% % 
% % clf; patch('Vertices',X1,'Faces',faces,'FaceVertexCData', tn,...
% %     'FaceColor','flat');axis square;
% % view(-70,4);colorbar;
% % zoom(1);drawnow;title('New object topology');
% % 
