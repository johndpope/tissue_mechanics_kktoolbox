function [R da dv dE X1 X2] = shape_distance(X1, X2)
%% calculate the shape distance between X1 and X2 
%% for the unnormalized version use object_distance (not written yet)

X1 = tr_inv(X1);
X2 = tr_inv(X2);

R = 1-(sum((abs(X1)-abs(X2)).^2)/(sum(X1.^2+X2.^2)));

[A1 V1 E1] = area_shps_notri(X1, 30);
[A2 V2 E2] = area_shps_notri(X2, 30);
da = abs(A1-A2);
dv = abs(V1-V2);
dE = abs(E1-E2);


disp([R da dv dE]);
figure;plot_sh_notri(X1,160);
figure;plot_sh_notri(X2,160);
