function plot_torus(zet,ita,phi)
%% given the toroidal coordinates plot the torus
% % %%% Convert from toroidal to explicit (parametric) Cartesian
% representation of torus
%load torus_frbf_3842;   % load X and F for a torus
load torus_frbf_24090;
a = 1;
x = a .* sinh(zet).* cos(phi)./(cosh(zet) - cos(ita));
y = a .* sinh(zet).* sin(phi)./(cosh(zet) - cos(ita));
z = a .* sin(ita)./(cosh(zet)-cos(ita)) .* ones(size(phi));
plot_mesh([x(:) y(:) z(:)],F);