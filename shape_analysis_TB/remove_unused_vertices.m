function [X2,F2] = remove_unused_vertices(X,F,R)
%%% remove the excess vertices listed in R, which is sorted
ind = 1:length(X); X2 = X;        %
for ix = 1:length(R),       % loop over the removed vertices
    r = R(1);              % this is the index to be removed
    X2 = X2(ind(1:end-ix+1)~=r,:);
    F(F>r) = F(F>r)-1;
    R = R(2:end)-1;%clip R by one.
end
F2 = F;