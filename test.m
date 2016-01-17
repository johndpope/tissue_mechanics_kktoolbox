n = 512; 
lev = 5; 
wav = 'db1';

% Wavelet based matrix multiplication by a vector: 
% a "good" example 
% Matrix is magic(512) Vector is (1:512)

m = magic(n); 
v = (1:n)'; 
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wav);

% ordinary matrix multiplication by a vector. 
p = m * v; 

% The number of floating point operations used is 524,288

% Compute matrix approximation at level 5. 
sm = m;
for i = 1:lev 
    sm = dyaddown(conv2(sm,Lo_D),'c'); 
    sm = dyaddown(conv2(sm,Lo_D'),'r'); 
end 

% The number of floating point operations used is 2,095,154

% The three steps: 
% 1. Compute vector approximation. 
% 2. Compute multiplication in wavelet domain. 
% 3. Reconstruct vector approximation.

sv = v; 
for i = 1:lev, sv = dyaddown(conv(sv,Lo_D)); end 
sp = sm * sv; 
for i = 1:lev, sp = conv(dyadup(sp),Lo_R); end 
sp = wkeep(sp,length(v));

% Now, the number of floating point operations used is 9058

% Relative square norm error in percent when using wavelets. 
rnrm = 100 * (norm(p-sp)/norm(p))