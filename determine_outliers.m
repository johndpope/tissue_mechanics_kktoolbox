function out = determine_outliers(c3, fac)
if nargin == 1,fac = 8;end
mu3 = mean(c3); % Data mean
sigma3 = std(c3); % Data standard deviation
out = find((c3 - mu3) > fac * sigma3);