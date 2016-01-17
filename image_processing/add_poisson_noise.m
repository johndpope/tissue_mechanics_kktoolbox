function I = add_poisson_noise(I,fac)
%% add some photon shot noise to image I
%%% add Poisson noise if needed
thresh = 0.1;
I(mat2gray(I)<thresh) = fac*max(I(:));;
I = poissrnd(I);