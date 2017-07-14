# tissue_mechanics_kktoolbox
Matlab code for simulating tissue mechanics
This is my sandbox for trying out new/alternative methods for simulating tissue mechanics. Please do not pull/fork/cite/rely on any of this code.

Example to generate an object and calculate its energy
-------------------------------------------------------
In Matlab
```json
s = shp_surface('discocyte');  % initializes an shp object, populates the coefficient vector with discocyte-coefficients and calculates the energy
figure;plot_pretty(s); % look at the surface
figure; plot_H(s); % plots the surface colored with local mean curvature
disp(s); % lists member variables of s including its bending energy
```
