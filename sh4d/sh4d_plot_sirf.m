function [f] = sh4d_plot_sirf(DNLK, dim, maxD)
% visualization of an intensity volume encoded in the coefficients D of DNLK
% which represent 4D spherical harmonics as in Matheny and Goldgof 1995
% i.e. 4D SH synthesis is done here
% Input: DNLK is a struct with fields, D, N, L, K and tau
% where D is a vector of coefficients
% N, L, K are vectors of equal size to D and correspond to a distinct basis
% function, tau is the period.
% dim: a vector of three values [xdim ydim zdim] with the size of the
% output intensity image I in pixels.
% Output: uint16 intensity image I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, dim = 32; end;
if nargin<3, maxD = 100;end;

if maxD>length(DNLK.D), maxD = length(DNLK.D);end
D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;
tau = DNLK.tau;
thresh = 1e-16;
av_i = round(dim/2);av_j = round(dim/2);av_k = round(dim/2);
%% generate theta, phi and r from the image
[ix iy iz] = ind2sub([dim dim dim],1:dim^3);
x = ix -av_i;y = iy-av_j;z = iz-av_k;  % move coordinate system to center of mass
[t p r] = kk_cart2sph(x,y,z);
f = zeros(dim^3,1,'single');



%%%% Let us truncate according to amplitude values D
[ds, IX] = sort(D.^2, 'descend');
D = D(IX(1:maxD)); 
N = N(IX(1:maxD));
L = L(IX(1:maxD));
K = K(IX(1:maxD));


%% generate the required basis vectors
s = length(L); %% number of functions in expansion
n = length(r(:)); %% number of data points
%%% Calculate the basis matrix using one processor
h = waitbar(0,'Generating basis matrix ....');
for S = 1:s,
    if abs(D(S))>thresh
        A = sh4d_basis_gen(N(S),L(S),K(S),p(:)',t(:)',r(:)', tau)';
        f = f+A*D(S);
    end
    waitbar(S/s);
end
close(h);

f = reshape(f,[dim dim dim]);
figure;kk_montage(mat2gray(f));impixelinfo;








