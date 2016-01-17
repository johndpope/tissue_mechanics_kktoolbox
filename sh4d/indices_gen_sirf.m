function [l,m] = indices_gen_sirf(L_max)
%prepare the indices for the SHP4D calculation  

l = [];
m = [];
for lix = 0:L_max,
    for kix = -lix:lix,
        l = [l lix];
        m = [m kix];
    end
end
