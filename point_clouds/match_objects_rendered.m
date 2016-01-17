function [r c all_Xdr cm m] = match_objects_rendered(all_X0, all_Xd, dangle)
%% match the objects defined by the pair of coordinate sets rotated with the specified
%% rotation (angle) in degrees

%[all_X0] = translate_XF(all_X0, -center_of_mass_set(all_X0));
C0 = center_of_mass(all_X0);

[all_Xdr, cm] = rotate_XF(all_Xd, [0 0 dangle]); % rotate the "all_Xd" set of objects to the angle of all_Xo around the center of mass
%[all_Xdr] = translate_XF(all_Xdr, -center_of_mass_set(all_Xdr));
Cdr = center_of_mass(all_Xdr);

%%%%%%%%%%%%%% Now let us overlap the two by minimizing the 
%%%%%%%%%%%%%% sum of squared minimum distances
global X Xdata
X = C0;
Xdata = Cdr;
Po = zeros(6,1);    % six parameters 3 translation and 3 rotation (Euler z-x-z convention)
                    % we start out with zeros because we should in fact be
                    % close to the transformation to all_X0
[P,resnorm,residual] = lsqnonlin(@register_lmobj,Po,[],[],...
    optimset('Algorithm', 'levenberg-marquardt',...
             'LargeScale','off', ...
             'MaxFunEvals', 1000,...
             'Display', 'on'));
[P,fval] = fminsearch(@register_obj, P, ...
    optimset('Display','off', 'MaxIter', 2000, 'TolX', 1e-3 ));
Cdr = translate_rotate_XF(Xdata, P);

%%% generate distance matrix between all points
% let's assume the fixed (model) points to be indexed along the direction
% of increasing row index
C0x = C0(:,1);C0x = C0x(:,ones(size(Cdr,1),1));
C0y = C0(:,2);C0y = C0y(:,ones(size(Cdr,1),1));
C0z = C0(:,3);C0z = C0z(:,ones(size(Cdr,1),1));

Cdrx = Cdr(:,1)';Cdrx = Cdrx(ones(size(C0,1),1), :);
Cdry = Cdr(:,2)';Cdry = Cdry(ones(size(C0,1),1), :);
Cdrz = Cdr(:,3)';Cdrz = Cdrz(ones(size(C0,1),1), :);

Rsq = ((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances squared
md = min(Rsq);  % row vector of  minimum distances squared

% let's match the objects based on the minimum distance
mdmx = md(ones(length(C0),1),:);
indx = find((mdmx(:)-Rsq(:))==0);
[r c] = ind2sub(size(Rsq), indx);   % r and c contain the matching object indices (into all_X0 and all_Xd respectively).
[r IX] = sort(r);                   % sort r (indices of all_X0) ascendingly
c = c(IX);                          % match the indices of all_Xd
d = md(c)';                          % match the distances to correspond to the r c pairs

%% make sure that we select only the object pairs that have a minimum
%% distance. This means that is the same object in one view is suggested to
%% match with more than one in another view, we select the one closest
for ix = 1:size(r,1),
   if sum(r==ix)>1,     % then we have multiple matches
       indx = find(r==ix);  % get the indices of the conflicting matches
       [B IX] = sort(d(indx));
       indx = indx(IX);
       keepr = r(indx(1));
       keepc = c(indx(1));
       r(indx) = [];c(indx) = [];
       r = [r; keepr];c = [c; keepc];
   end
end
[r IX] = sort(r);c = c(IX);d = md(c)';  

[r, m, n] = unique(r,'rows','first'); c = c(m);% r and c contain the matching object indices (into all_X0 and all_Xd respectively).
% d = md(c)';
































