function [X Q F] = E_XQ_gen()
X(1,:) = [-5.5 8 1];
X(2,:) = [-1.5 8 1];
X(3,:) = [5.5 8 1];
X(4,:) = [-1.5 5.5 1];
X(5,:) = [5.5 5.5 1];
X(6,:) = [-1.5 1.5 1];
X(7,:) = [-1.5 -1.5 1];
X(8,:) = [-1.5 -5.5 1];
X(9,:) = [-1.5 -8 1];
X(10,:) = [5 1.5 1];
X(11,:) = [5 -1.5 1];
X(12,:) = [5.5 -5.5 1];
X(13,:) = [5.5 -8 1];
X(14,:) = [-5.5 -8 1];
X2 = X;X2(:,3) = X2(:,3)-2;
X = [X;X2];
F(1,:) = [1 2 9 14];
F(2,:) = [2 3 5 4];
F(3,:) = [6 10 11 7];
F(4,:) = [8 12 13 9];
F(5:8,:) = F(1:4,:)+14; % generate the bottom faces
F(end+1,:) = [14 1 15 28];
F(end+1,:) = [14 13 27 28];
F(end+1,:) = [1 15 17 3];
 F(end + 1,:) = [3 5 19 17];
F(end+1,:) = [10 24 25 11];
F(end+1,:) = [13 12 26 27];
F(end+1,:) = [5 4 18 19];
F(end+1,:) = [4 7 21 18];
F(end+1,:) = [6 10 24 20];
F(end+1,:) = [11 7 21 25];
F(end+1,:) = [7 8 22 21];
F(end+1,:) = [8 12 26 22];
Q = F;
if nargout == 3
    F = quad2tri(Q);
end
function F = quad2tri(Q)
%% generate F by taking a diagonal through every quad
F = [];
for(qix = 1:size(Q,1))
    q = Q(qix,:);
    f1 = [q(1) q(2) q(3)];
    f2 = [q(3) q(4) q(1)];
    F = [F;f1;f2];
end
% 
% %% plot
% clf;edge_color = 0.0;
% patch('vertices',X,'faces',F,'facecolor','green','edgecolor',[edge_color edge_color edge_color]);
% axis on;axis equal;cameramenu;