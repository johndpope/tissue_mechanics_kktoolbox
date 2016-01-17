function [all_X_new] = translate_XF(all_X, cm)
%% translate the group of shapes 
convert_to_cell = 0;
if ~iscell(all_X), convert_to_cell = 1;X{1} = all_X; all_X = X;end
X = cell2mat(all_X');   % concatenate all X
X = [X ones(size(X,1), 1)];
% generate the translation matrix
T = [1 0 0 cm(1);...
    0 1 0 cm(2);...
    0 0 1 cm(3);...
    0 0 0  1];
% apply the transform
Xp = (T*X')';
% generate the new set of objects
counter = 1;
all_X_new = all_X;
for ix = 1:size(all_X,2),
    l = size(all_X{ix},1);
    all_X_new{ix} = Xp(counter:l+counter-1, 1:3);
    counter = counter + l;
end

if convert_to_cell, all_X_new = cell2mat(all_X_new);end