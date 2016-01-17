%%% generate a quadmesh two voxel and test the matrices as in C.B. 1995
close all;clc;
% % two-voxel C.B.
ix = 0;D.X(ix+1, :) = [ 9 14 6];D.N{ix+1} = [1  7  6  9  3 4]+1;
ix = 1;D.X(ix+1, :) = [ 10 14 6];D.N{ix+1} = [0  3  4  5  2  8  7  6]+1;
ix = 2;D.X(ix+1, :) = [ 11 14 6];D.N{ix+1} = [1  4  5  11  8  7]+1;
ix = 3;D.X(ix+1, :) = [9 15 6];D.N{ix+1} = [4  1  0  6  9  10]+1;
ix = 4;D.X(ix+1, :) = [10 15 6];D.N{ix+1} = [3  9  10  11  5  2  1  0]+1;
ix = 5;D.X(ix+1, :) = [11 15 6];D.N{ix+1} = [4  10  11  8  2  1]+1;
ix = 6;D.X(ix+1, :) = [ 9 14 7];D.N{ix+1} = [7  10  9  3  0  1]+1;
ix = 7;D.X(ix+1, :) = [ 10 14 7];D.N{ix+1} = [6  0  1  2  8  11  10  9]+1;
ix = 8;D.X(ix+1, :) = [ 11 14 7];D.N{ix+1} = [7  1  2  5  11  10]+1;
ix = 9;D.X(ix+1, :) = [ 9 15 7];D.N{ix+1} = [10  4  3  0  6  7]+1;
ix = 10;D.X(ix+1, :) = [10 15 7];D.N{ix+1}= [9  6  7  8  11  5  4 3]+1;
ix = 11;D.X(ix+1, :) = [11 15 7];D.N{ix+1} = [10  7  8  2  5  4]+1;
F = [0 1 7 6  ;1 2 8 7 ;7 8 11 10;6 7 10 9;2 5 11 8;...
     4 5 11 10;3 4 10 9;0 3 9 6  ;0 1 4 3;1 2 5 4];
F = F + 1;
X = D.X;
%F = [F(:,[1 2 3]); F(:,[1 4 3])];% convert to triangulation??
%[X,F] = subdivide(X,F,2);
ixN = 1;
ixS = 12;
[t,p, E, L, face_memb, A,b,Aq, bq, dl] = qvbijective_map_gen(X,F, ixN, ixS);