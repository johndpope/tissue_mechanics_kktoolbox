function [I, A] = sh4d_plot(DNLK, dim, A, flag)
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
if nargin <3, A = []; flag = 1;end
if nargin <4, flag = 1;end
D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;
tau = DNLK.tau;
thresh = 1e-20;

I = zeros([dim dim dim],'single');

av_i = round(size(I,1)/2);av_j = round(size(I,2)/2);av_k = round(size(I,3)/2);

%% generate theta, phi and r from the image
[ix iy iz] = ind2sub(size(I),1:length(I(:)));
x = ix -av_i;y = iy-av_j;z = iz-av_k;  % move coordinate system to center of mass
[t p r] = kk_cart2sph(x,y,z);

%% generate the required basis vectors
% if isempty(A),

    m = length(L); %% number of functions in expansion
%     disp(m);
    n = length(r(:)); %% number of data points
    %%% Calculate the basis matrix using one processor
    A  = zeros(n, m, 'single');
    h = waitbar(0,'Generating basis matrix ....');
    for S = 1:m,
        if abs(D(S))>thresh,
        A(:, S) = sh4d_basis_gen(N(S),L(S),K(S),p(:)',t(:)',r(:)', tau)';
        else
            A(:, S) = 0;
        end
        waitbar(S/m);
    end
    close(h);
% end

f = A*D;
indx = sub2ind(size(I),ix, iy, iz);

I(indx) = f;
if flag, figure;kk_montage(mat2gray(I));impixelinfo;end








