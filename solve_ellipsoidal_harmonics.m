function [xp yp xc yc A] = solve_ellipsoidal_harmonics(Jmax, X)
%% solve the 2D outline problem, given that X is ordered squentially
%%% if X is not ordered sequentially then we need generate a sorted/ordered
%%% X for instance using the  code (or similar) at the end

t = 0:2*pi/length(X):(2*pi-2*pi/length(X));
A = ellipsoidal_basis_gen(Jmax, t);
% construct the data vectors
bx = X(:,1);
by = X(:,2);
% solve system of equations
xc = bx\A;
yc = by\A;
%
xp = A*xc';
yp = A*yc';


%%
% % % % % %% order the points according to their neighborhood relationship
% % % % % close all;
% % % % % X = zeros(length(x),2);
% % % % % ixt = zeros(length(x),1);
% % % % % counter = 1;
% % % % % ixc = 1561;
% % % % % ixt(ixc) = 1;
% % % % % X(counter,:) = [x(ixc) y(ixc)];
% % % % % while counter<=length(x)
% % % % %     disp(counter);
% % % % %     mind = inf;
% % % % %     for jx = 1:length(x),
% % % % %         d = sqrt(sum(([x(ixc) y(ixc)]-[x(jx) y(jx)]).^2));
% % % % %         %if d~=0
% % % % %             if ixt(jx)==0,
% % % % %                 if d<mind
% % % % %                     mind = d;
% % % % %                     ixn = jx;   % set the next point
% % % % %                     X(counter,:) = [x(jx) y(jx)];
% % % % %                 end
% % % % %             end
% % % % %         %end
% % % % %     end
% % % % %     %plot([x(ixc) x(ixn)], [y(ixc) y(ixn)], '-');hold on; drawnow;
% % % % %     ixc = ixn;      % set the new current point
% % % % %     ixt(ixc) = 1;   % mark point as taken
% % % % %     counter = counter + 1;
% % % % % end
% % % % % X(end,:) = [179 323];
% % % % % tmp = X(:,1);
% % % % % X(:,1) = X(:,2);
% % % % % X(:,2) = tmp;
% % % % % % transfer to center of mass
% % % % % cm = [sum(X(:,1)) sum(X(:,2))]/length(X);
% % % % % X(:,1) = X(:,1)-cm(1);
% % % % % X(:,2) = X(:,2)-cm(2);
% % % % % save X X;