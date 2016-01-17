function ix = get_LK_index(L,K)
% assuming our usual arrangement of indices L and K, we get the linear
% index reflecting the position of the L,K pair

[l,m, flags] = indices_gen(ones((max(L)+1)^2, 1));

lcol = (l==L);
mcol = (m==K);
ix = find(sum([lcol(:) mcol(:)],2)==2);