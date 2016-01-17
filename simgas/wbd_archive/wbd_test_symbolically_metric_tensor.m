%%% symbolically obtain right Cauchy green tensor using 2 different
%%% approaches and show that they are equivalent
clc;clear all;
syms xt yt zt xp yp zp real
syms Xt Yt Zt Xp Yp Zp real

m1 = [xt yt zt];
m2 = [xp yp zp];
m3 = cross(m1,m2)/sqrt(sum(cross(m1,m2)));
g11 = dot(m1,m1);
g12 = dot(m1,m2);
g21 = g12;
g22 = dot(m2,m2);
gij =[g11 g21;g12 g22]; 
mu1 = cross(m2,m3)/sqrt(sum(cross(m1,m2)));% pretty(simplify(mu1(1)));
mu2 = cross(m3,m1)/sqrt(sum(cross(m1,m2)));

%% now for the deformed
M1 = [Xt Yt Zt];
M2 = [xp yp zp];
M3 = cross(M1,M2)/sqrt(sum(cross(M1,M2)));
G11 = dot(M1,M1);
G12 = dot(M1,M2);
G21 = G12;
G22 = dot(M2,M2);
Gij =[G11 G21;G12 G22]; 
Mu1 = cross(M2,M3)/sqrt(sum(cross(M1,M2)));
Mu2 = cross(M3,M1)/sqrt(sum(cross(M1,M2)));
%% right Cauchy green tensor and strain tensor for shear only as in AMoS
C1 = kron(M1',mu1) + kron(M2',mu2);
eps1 = 1/2*((C1'*C1)-eye(size(C1)));
pretty(simplify(eps1));
%% strain tensor method II directly from metric tensors
gamma11 = 1/2*(G11-g11)./sqrt(g11)./sqrt(g11);
gamma12 = 1/2*(G12-g12)./sqrt(g11)./sqrt(g22);
gamma21 = 1/2*(G12-g12)./sqrt(g11)./sqrt(g22);
gamma22 = 1/2*(G22-g22)./sqrt(g22)./sqrt(g22);
eps2 = [gamma11 gamma21 0;gamma12 gamma22 0; 0 0 0];
C2 = 2*eps2+eye(size(eps2));

disp(eig(C1)-eig(C2));














