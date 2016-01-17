function [R] = shape_distance(X1, X2)
%% calculate the shape distance between X1 and X2 
%% for the unnormalized version use object_distance (not written yet)
%% at the moment the shape destance is based on the first few coefficients only

L_max = get_L_max(X1);

X1 = s_inv(tr(X1,L_max));
X2 = s_inv(tr(X2,L_max));
[xc yc zc] = get_xyz_clks(X1);
[xc1 yc1 zc1] = get_xyz_clks(X2);


% R = sum((X1(:)-X2(:)).^2);

vec = 2:length(xc);
Ar = [xc(vec) yc(vec) zc(vec)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD1 = sqrt(diag(D));

Ar = [xc1(vec) yc1(vec) zc1(vec)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD2 = sqrt(diag(D));
R = sum((sqrtdD1-sqrtdD2).^2);



% % % L_max = 3;
% % % %%% for L = 1
% % % vec = 2:4;
% % % Ar = [xc(vec) yc(vec) zc(vec)];
% % % [V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
% % % sqrtdD1 = sqrt(diag(D));
% % % 
% % % Ar = [xc1(vec) yc1(vec) zc1(vec)];
% % % [V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
% % % sqrtdD2 = sqrt(diag(D));
% % % R = sum((sqrtdD1-sqrtdD2).^2);
% % % 
% % % %%% for L = 2
% % % A = [xc(5:9) yc(5:9) zc(5:9)];
% % % [V D] = eig(A'*A);  %% diagonalize L = 2
% % % D21 = sqrt(diag(D));
% % % A = [xc1(5:9) yc1(5:9) zc1(5:9)];
% % % [V D] = eig(A'*A);  %% diagonalize L = 2
% % % D22 = sqrt(diag(D));
% % % R = R + sum((D21-D22).^2);
% % % 
% % % %%% for L = 3
% % % vec = 10:16;
% % % A = [xc(vec) yc(vec) zc(vec)];
% % % [V D] = eig(A'*A);  %% diagonalize L = 2
% % % D21 = sqrt(diag(D));
% % % A = [xc1(vec) yc1(vec) zc1(vec)];
% % % [V D] = eig(A'*A);  %% diagonalize L = 2
% % % D22 = sqrt(diag(D));
% % % R = R + sum((D21-D22).^2);
