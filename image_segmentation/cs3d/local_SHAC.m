function [segmented mosaic] = local_SHAC(I, mosaic, num, lbnd)
%%% loop through the mosaic region, pick in each region the maximum
%%% intensity, fit an SHP local active contour to the region around this
%%% maximum intensity.