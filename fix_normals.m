function [Xout Fout] = fix_normals(X,F)
fvin.vertices = X;fvin.faces = F;
[coordNORMALS,fvout] = COMPUTE_mesh_normals(fvin);%% fix the normal directions
Xout = fvout.vertices;
Fout = fvout.faces;