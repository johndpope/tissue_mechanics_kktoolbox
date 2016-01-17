%%%%%% load a torus mesh and convert the Cartesian coordinates to toroidal
%% note: Toroidal coordinates notation as in Milligan and Fraguas 1994
clc;close all;
%load torus_frbf_3842;   % load X and F for a torus
load torus_frbf_24090;
x = X(:,1);y = X(:,2);z = X(:,3);
[phi,R,Z] = cart2pol(x, y, z);   % convert mesh to cylindrical coordinates angle theta, radius r and height z

% convert from cylindrical to toroidal
zet = 2;
ita = atan2(Z.*sinh(zet),(R.*cosh(zet)- sinh(zet)));
plot_torus(zet,ita,phi)
%% %%%%%%%%% Calculate the Toroidal Harmonics
N = 10;
M = 10;
A = zeros(M,N);
B = zeros(M,N);
a = zeros(N,1);
b = zeros(N,1);
c = zeros(M,1);
d = zeros(M,1);

A(1) = 0;
B(10,10) = 1e-1;
a(1) = 0;
b(10) = 1e-1;
c(10) = 0.01;
d(10) = 0.01;

S = 0;
P0 = P_n_m(0,0,zet);
Q0 = Q_n_m(0,0,zet);
fac = sqrt(cosh(zet)-cos(ita));
for mix = 0:M-1,
    for nix = 0:N-1,
        m = mix+1;
        n = nix+1;
        t1 = (A(m,n).*P_n_m(n,m,zet)/P0 + B(m,n).*Q_n_m(n,m,zet)/Q0);
        t2 = (a(n)*cos(n*ita) + b(n)*sin(n*ita)).*(c(m)*cos(m*phi) + d(m)*sin(m*phi));
        disp([m n t1 sum(t2)]);
        S = S + fac.*(t1.*t2);
    end
end
disp(sum(abs(S)));
% % plot the surface as a color coding on a torus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = [x y z];
% dfig(1);clf;colormap hot;patch('Vertices', X, 'Faces', F,'FaceVertexCData',S, 'FaceColor', 'interp', 'EdgeColor','none');
% daspect([1 1 1]);axis on;lighting gouraud;
% view(2);rotate3d; drawnow;xlabel('x');ylabel('y');zlabel('z');camlight

%plot_mesh([x(:) y(:) z(:)],F);
