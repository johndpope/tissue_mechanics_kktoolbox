function Ir = image_rotate(I,angd, flag)
%%% rotates the 3D image around an axis defined by flag 1 = z, 2 = y, 3=x
%%% by an angle in degrees!! angd.

%%% to rotate around anything else but x axes (SPIM and MATLAB both define
%%% the same x axis) permutations have to be done first


%%%
Ir = [];
if angd==0,
    Ir = I;
else
    for ix = 1:size(I,3),
       Ir(:,:,ix) = imrotate(I(:,:,ix), deg2rad(angd));
    end 
end