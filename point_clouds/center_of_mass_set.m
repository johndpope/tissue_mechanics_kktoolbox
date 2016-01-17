function cm = center_of_mass_set(all_X)
%% return the 3-vector for the center of mass all objects in all_X
X = cell2mat(all_X');   % concatenate all X
X = [X ones(size(X,1), 1)];
cm = sum(X,1)./size(X,1);