function [t,p, shear, area_res] = uniform_dist_gen(t,p,S, F,L,E, F_areas_relative, maxiter, DMflag, JacPatCalc, AreaCalcEqual)
%% Calculates the uniform configuration on the sphere, based
%% on the triangulation F and L, and starting at t and p
disp('Optimizing the sampling on the sphere');
if nargin <8, maxiter=50; DMflag = 0;JacPatCalc = 1;AreaCalcEqual = 0;end
% niter = 2000;az = -95;el = 16;nvert = length(t);
%[t,p] = relax_mesh(t, p, E, F,1000, az, el); %%%%% let the mesh relax    
[t,p, shear, area_res] = optimize_set(F, S, t, p, F_areas_relative,maxiter, DMflag, JacPatCalc, AreaCalcEqual);    % do an optimization

%%% now that we have the highly sheared triangles we can modify the
%%% triangulation of the shape and remap.