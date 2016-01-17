function level = graythresh_sct(I)
% compute a threshold level using the stable count thresholding algorithm
% of Russel et al. Biophysical J. 2009.
% This method is supposed to be superior to Otsu's method
% Input: I: 3D image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        delta: used for central difference calculation to estimate asymptote: set to 10
%        ad   : factor in denominator for calculation of alpha to decide
%        when the smallest threshold is reached
thresh = 1e1;
I = mat2gray(I);
T_vec = [];
NT_vec = [];
for T = 0:0.01:1
    NT      = sum(I(:)>T);  
    T_vec = [T_vec T];
    NT_vec = [NT_vec NT];
end

indx = find((abs(diff(NT_vec))<thresh));
level = T_vec(indx(1));









