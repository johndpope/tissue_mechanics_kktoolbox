function [X1, faces] = kk_subdivide(X, F)
%% subdivide the triangles in F by taking the center of mass of the triangle
%% to be the new vertex. The result of the subdivision is that every triangle gives three new ones
faces = [];
X1 = X;
for ix = 1:size(F,1),
    newV = [(X(F(ix,1),1)+  X(F(ix,2),1) + X(F(ix,3),1))/3 ...
            (X(F(ix,1),2) + X(F(ix,2),2) + X(F(ix,3),2))/3 ...
            (X(F(ix,1),3) + X(F(ix,2),3) + X(F(ix,3),3))/3];
        X1 = [X1;newV];
        newix = length(X1);
        f1 = [F(ix,1) F(ix,2) newix];
        f2 = [F(ix,2) F(ix,3) newix];
        f3 = [F(ix,3) F(ix,1) newix];
        faces = [faces;f1;f2;f3];
end