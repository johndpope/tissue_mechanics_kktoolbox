function [d] = shp_d(X1, X2)
%% calculate the shape distance between two r_inv shp shapes
d = 1/sqrt(4*pi) * sqrt(sum((X1-X2).^2));
