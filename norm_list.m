function n = norm_list(V)
% returns a one column list of norms for the vectors (rows) given in V
n = zeros(size(V, 1), 1);
for ix = 1:size(V,1)
    n(ix)= norm(V(ix,:));
end
