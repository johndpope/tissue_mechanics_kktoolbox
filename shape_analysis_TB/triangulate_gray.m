function [X, F, L,E] = triangulate_gray(I,x_pixel,y_pixel,z_pixel,thres,nplanes,plotflag)
%% Triangulate the shape in I based on the isosurface method in matlab
%% according to the threshold thres
X = []; F = []; L = []; E = [];
verbose = 1;
I2 = I;
medf = 0;
if medf, I2 = smooth3(I,'box',[medf medf medf]);end    % optional
[F,X] = isosurface(I2,thres); % optimal threshold still needs to be det.
%disp('Isosurface calculation and triangulation');
%save data_intermediate_triangulategray
% for ix = 1:rflag
% [X F] =subdivide(X,F,1); %%% refine mesh if necessary
% end

%disp(['Number of vertices = ' num2str(length(X)) '  Number of faces = ' num2str(length(F))]);
%if verbose, patch('Vertices', X, 'Faces', F, 'FaceColor', 'blue'); axis equal;end
if  length(X)*2-4~=length(F),disp([length(X)*2-4 length(F) (length(X)*2-4)-length(F)]);
    X = [];
    disp('Shape is not closed');
end

if nargout >2
    nvert = length(X);
    E = zeros(nvert,2);counter = 0;
    % determine the links L
    for ix = 1:length(X),   % loop over the vertices
        if ~mod(ix,2000),if plotflag,disp(ix);end;end; fmemb = ismember(F, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
        [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
        links = [];
        for ik = 1:length(ig),
            links = [links F(ig(ik),:)];
        end
        L{ix} = unique(links(links~=ix));   % only record the links that are not ix
        % create the list of all edges
        llinks = L{ix};
        for jx = 1:length(llinks),
            counter = counter + 1;
            E(counter,:) = [ix llinks(jx)];
        end;
    end;
end
% % center the object at the center of mass
% % if ~isempty(X)
% %     x = X(:,1) - mean(X(:,1)); y = X(:,2)-mean(X(:,2)); z = X(:,3)-mean(X(:,3));
% %     x = x*x_pixel; y = y*y_pixel;z = z*z_pixel;X = [x y z];
% % end
%% Look at the shape
if plotflag && verbose
    figure
    %axis square;patch('Vertices',X,'Faces',F,'FaceVertexCData',hsv(size(F,1)),'FaceColor','flat', 'FaceAlpha', .8);
    axis square;axis equal; axis tight;xlabel('microns'); ylabel('microns');zlabel('microns');
    [xs,ys,zs,D] = subvolume(I,[nan,nan,nan,nan,1,nplanes]);
    xs = xs*x_pixel; ys = ys*y_pixel;zs = zs*z_pixel;
    p1 = patch(isosurface(xs,ys,zs,D, thres),...
        'FaceColor','red','EdgeColor','none');
    isonormals(xs,ys,zs,D,p1);
    p2 = patch(isocaps(xs,ys,zs,D, thres),...
        'FaceColor','interp');%,'EdgeColor','none');
    view(3); axis tight; daspect([1,1,1])
    axis equal; axis tight;axis off
    colormap(gray(100))
    camlight right; camlight left; lighting gouraud;lightangle(100,-90)
    view(3);drawnow;rotate3d; hold off
end
