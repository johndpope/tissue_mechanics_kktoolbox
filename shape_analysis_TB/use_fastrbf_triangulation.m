% plotflag = 1;   
nvert = length(X);
plotflag = 0;
    E = zeros(nvert,2);counter = 0;
    L = {};
    % determine the links L
    for ix = 1:length(X),   % loop over the vertices
        if ~mod(ix,1000),if plotflag,disp([num2str(ix) ' of ' num2str(length(X))]);end;end; fmemb = ismember(F, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
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
    % % center the object at the center of mass
    %%save data_intermediate_mesh_info
x = X(:,1) - mean(X(:,1)); y = X(:,2)-mean(X(:,2)); z = X(:,3)-mean(X(:,3));
%x = x*x_pixel; y = y*y_pixel;z = z*z_pixel;X = [x y z];
if (length(X)-length(E)/2+length(F)~=2),
    disp((length(X)-length(E)/2+length(F)));
    euler_violation = 1;
%     error('Euler characteristic violation');
else
    euler_violation = 0;
end
