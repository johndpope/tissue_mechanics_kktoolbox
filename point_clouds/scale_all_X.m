function [all_X, cm] = scale_all_X(all_X, scale)
%% scale the group of shapes around each center of mass individually
%% all_X is a cell array

for gix = 1:length(all_X)
    X = all_X{gix};
%% center of mass
cm = sum(X,1)./length(X);cm = cm(:)';
% translate the center of mass to the origin
X = X-cm(ones(1,size(X,1)), :);
% apply the scaling
Xp = (X*scale);
% translate back the center of mass to where it was
Xp = Xp+cm(ones(1,size(X,1)), :);
all_X{gix} = Xp;

end
