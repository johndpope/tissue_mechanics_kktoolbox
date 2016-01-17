function T = generate_T(X,F)
%% generate properties per vertex
[Edge L face_memb] = edge_info_02(X,F);
[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA,normals] = triangulated_props(X, F, 0);
T.Hv = H;
T.Kv = K;
T.Nv = N;
x = X(:,1);y = X(:,2);z = X(:,3);
for cix = 1:length(X)        % loop over indices
    f = face_memb{cix};
    acum = 0;
    dcum = 0;
    for(fix = 1:length(f))  % loop over those faces
        V1 = [x(F(f(fix),1)) y(F(f(fix),1)) z(F(f(fix),1)) ];
        V2 = [x(F(f(fix),2)) y(F(f(fix),2)) z(F(f(fix),2)) ];
        V3 = [x(F(f(fix),3)) y(F(f(fix),3)) z(F(f(fix),3)) ];
        [a n d] = tri_prop(V1, V2, V3);
        acum = acum + a;
        dcum = dcum + (1-max(d)/min(d));
    end
    acum = acum/numel(f);
    dcum = dcum/numel(f);
    T.Av(cix) = acum;
    T.dv(cix) = dcum;
end
T.Av = T.Av(:)/sum(T.Av(:));
T.dv = T.dv(:);
%% determine properties per face
T.Af = F_areas/sum(F_areas);
for(fix = 1:size(F,1))
    T.df(fix) = (T.dv(F(fix,1)) + T.dv(F(fix,2)) +T.dv(F(fix,3)))/3;
    T.Hf(fix) = (T.Hv(F(fix,1)) + T.Hv(F(fix,2)) +T.Hv(F(fix,3)))/3;
    T.Kf(fix) = (T.Kv(F(fix,1)) + T.Kv(F(fix,2)) +T.Kv(F(fix,3)))/3;
end
T.df = T.df(:);
T.Hf = T.Hf(:);
T.Kf = T.Kf(:);
T.Nf = normals;
