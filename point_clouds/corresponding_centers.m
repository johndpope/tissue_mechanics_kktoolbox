function [X Xr] = corresponding_centers(X,Xr,n)
%%%determine the corresponding centers using local geometric descriptors

%%%%%% calculate the distance matrix for the first dataset and construct
%%%%%% neighbor list
%%% n is the number of neighboring beads to consider when building the local descriptor
        C0x = X(:,1);C0x = C0x(:,ones(size(X,1),1));
        C0y = X(:,2);C0y = C0y(:,ones(size(X,1),1));
        C0z = X(:,3);C0z = C0z(:,ones(size(X,1),1));

        Cdrx = X(:,1)';Cdrx = Cdrx(ones(size(X,1),1), :);
        Cdry = X(:,2)';Cdry = Cdry(ones(size(X,1),1), :);
        Cdrz = X(:,3)';Cdrz = Cdrz(ones(size(X,1),1), :);

        Rd = sqrt((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances
        nlist1 = {}; % neighbor list of length "length(X)"
        for row = 1:size(Rd,1),      % loop over the beads of image 1
            bead = X(row,:);    % get the bead coordinates
            d = Rd(row,:);      % get the row vector of distances
            %d(row,row) = [];    % delete the distance with itself (zero)
            [d IX] = sort(d);       % sort d in ascending order and obtain the indices of the partners
            d = d(1:n+1);        % take the n first elements of d
            IX = IX(1:n+1);      % indices of the n closes neighbors (plus self)
            neigh = X(IX(2:end),:);    % list of coordinates of the n neighbors
            neigh = neigh-bead(ones(size(neigh,1),1),:);        % translate the point cloud to be centered around the bead
            nlist1{row} = neigh;
        end
        
%%%%%% calculate the distance matrix for the second dataset
        C0x = Xr(:,1);C0x = C0x(:,ones(size(Xr,1),1));
        C0y = Xr(:,2);C0y = C0y(:,ones(size(Xr,1),1));
        C0z = Xr(:,3);C0z = C0z(:,ones(size(Xr,1),1));

        Cdrx = Xr(:,1)';Cdrx = Cdrx(ones(size(Xr,1),1), :);
        Cdry = Xr(:,2)';Cdry = Cdry(ones(size(Xr,1),1), :);
        Cdrz = Xr(:,3)';Cdrz = Cdrz(ones(size(Xr,1),1), :);

        Rd = sqrt((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances
        nlist2 = {}; % neighbor list of length "length(X)"
        for row = 1:size(Rd,1),      % loop over the beads of image 1
            bead = Xr(row,:);    % get the bead coordinates
            d = Rd(row,:);      % get the row vector of disctances
            %d(row,row) = [];    % delete the distance with itself (zero)
            [d IX] = sort(d);       % sort d in ascending order and obtain the indices of the partners
            d = d(1:n+1);        % take the n first elements of d
            IX = IX(1:n+1);
            neigh = Xr(IX(2:end),:);    % list of coordinates of the n neighbors
            neigh = neigh-bead(ones(size(neigh,1),1),:);        % translate the point cloud to be centered around the bead
            nlist2{row} = neigh;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determination of the correspondence
gld = zeros(length(nlist1)*length(nlist2),3,'double');  % store the indices of the pair and the dissimilarity measure
counter = 0;
for ix = 1:length(nlist1),      % loop over list 1
    X1 = nlist1{ix};            % obtain the list of coordinates of the neighbors of ixth bead of image 1 relative to the bead
    disp([num2str(ix) ' of ' num2str(length(nlist1)) ' iterations']);
    for jx = 1:length(nlist2),  % loop over list 2
        counter = counter + 1;
        X2 = nlist2{jx};
        %% calculate the similarity using procrustes analysis
        [d,Z,tr] = procrustes(X1,X2);
        
        trUnscaled.T = tr.T;
        trUnscaled.b = 1;
        trUnscaled.c = mean(X1) - mean(X2) * trUnscaled.T;
        ZUnscaled = X2 * trUnscaled.T + repmat(trUnscaled.c,n,1);
        dUnscaled = sum((ZUnscaled(:)-X1(:)).^2) ...
            / sum(sum((X1 - repmat(mean(X1,1),n,1)).^2, 1));


%         gld(counter,:)= [ix jx dissimilarity(X1, ZUnscaled)];
        gld(counter,:)= [ix jx dUnscaled];
        tr.Z = Z;
        all_transforms{counter} = trUnscaled;
    end
end
[B,IX] = sort(gld(:,3), 'ascend');
gld = gld(IX,:);
all_transforms = all_transforms(IX);
nbeads = min(size(X,1), size(Xr,1));
X = X(gld(1:nbeads,1),:);
Xr= Xr(gld(1:nbeads,2),:);
all_transforms = all_transforms(1:nbeads);
disp(' ');
hist(gld(:,3),1000);
%%% extract the unique combinations such that each bead is used only once
%%% in combination with its best match
rm = [];
for ix = 1:length(Xr)
    vec = find(Xr(:,1)==Xr(ix,1));    % get row indices into gld which have multiple occurrences of vec.
    while length(vec)>1,
       rm = [rm;vec(end)];
       vec(end) = [];
    end
end
X(rm,:) = [];
Xr(rm,:) = [];
all_transforms(rm) = [];
%%%% look at the  transforms and identify false correspondences that have
%%%% outlier transforms
T = all_transforms{1}.T;
cutoff = 4.0;
rm = [];
for ix = 2:length(all_transforms)
    res = all_transforms{ix}.T-T;
    %disp(sum(abs(res(:))));
    if sum(abs(res(:)))>cutoff, rm = [rm;ix];end
end
X(rm,:) = [];
Xr(rm,:) = [];
all_transforms(rm) = [];

% % % uncomment to debug stuff
% % distance = sqrt((X(:,1)-Xr(:,1)).^2 + (X(:,2)-Xr(:,2)).^2 + (X(:,3)-X(:,3)).^2);
% % K = size(X,1);
% % figure;
% % figure; kk_plot3(X(1:K,:),'b');hold on;kk_plot3(Xr(1:K,:)+0.2, 'r');hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function d = dissimilarity(X1, X2)
% returns a measure of dissimilarity based on the closest point criterion
C0x = X1(:,1);C0x = C0x(:,ones(size(X2,1),1));
C0y = X1(:,2);C0y = C0y(:,ones(size(X2,1),1));
C0z = X1(:,3);C0z = C0z(:,ones(size(X2,1),1));

Cdrx = X2(:,1)';Cdrx = Cdrx(ones(size(X1,1),1), :);
Cdry = X2(:,2)';Cdry = Cdry(ones(size(X1,1),1), :);
Cdrz = X2(:,3)';Cdrz = Cdrz(ones(size(X1,1),1), :);

Rsq = ((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances squared
md = min(Rsq);  % row vector of  minimum distances squared
d  = sum(md(:));

































