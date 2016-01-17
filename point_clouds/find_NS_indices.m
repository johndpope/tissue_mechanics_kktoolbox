function [indxN indxS] = find_NS_indices(X,NP, SP)
indxN = 0;indxS = 0;
dNmin = inf;dSmin = inf;
for (ix = 1:length(X))
    v = X(ix,:);
    dN = sqrt(sum((v-NP).^2));
    dS = sqrt(sum((v-SP).^2));
    if dN<dNmin, dNmin = dN;indxN = ix;end
    if dS<dSmin, dSmin = dS;indxS = ix;end
end