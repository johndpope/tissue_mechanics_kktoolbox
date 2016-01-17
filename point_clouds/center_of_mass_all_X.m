function cm_vec = center_of_mass_all_X(all_X)
%% return the vector for the centers of mass of the objects in all_X
cm_vec = [];
for ix = 1:length(all_X),
    X = all_X{ix};
    cm_vec = [cm_vec;sum(X,1)./length(X)];
end