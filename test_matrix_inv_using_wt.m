%%% the idea is to approximate a dense (design) matrix A with its wavelet transform
%%% Aw, and use Aw for solving the system Aw . xw = bw , of course using a
%%% wavelet transformed data vector b. The solution x is found by inverse wavelet transformation.
clc
%% let us first make A, x and b
dim = 512;
A = rand(dim);
 %A  = mat2gray(imread('testmatrix.tif'));dim = size(A,1);
x = 1:dim; x = x(:);
b = A*x;

% disp(norm(x-A\b));
% Relative square norm error in percent when using wavelets. 
rnrm = 100 * (norm(x-A\b)/norm(x));
disp(rnrm);

%% now let us transform A and b
%% compute matrix approximation at level 5
lev = 3;%log2(dim)-1
wav = 'db5';
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wav);
sm = A;
for i = 1:lev 
    sm = dyaddown(conv2(sm,Lo_D),'c'); 
    sm = dyaddown(conv2(sm,Lo_D'),'r'); 
end

%disp(sm);

% The three steps: 
% 1. Compute vector approximation. 
% 2. Compute vector x in wavelet domain using QR decomposition. 
% 3. Reconstruct vector x approximation.

sb = b; 
for i = 1:lev, sb = dyaddown(conv(sb,Lo_D)); end 
sx = sm\sb; 
for i = 1:lev, sx = conv(dyadup(sx),Lo_R); end 
sx = wkeep(sx,length(b)); 
% Relative square norm error in percent when using wavelets. 
rnrm = 100 * (norm(x-sx)/norm(x));
disp(rnrm);

















