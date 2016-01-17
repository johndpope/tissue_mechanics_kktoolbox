function [clks, U, S, V, invS] = sh_expand(R, L_max, t,p, recalc)
%% Calculates the spherical harmonic analysis (expansion) of the function R
%% defined on the sphere at the points t, p. The expansion is calculated to
%% order L = L_max.
%% Note that R can also be a coordinate mapped onto the sphere
%% This is particularly important in the case of multiple
%% expansion.
%% USAGE: [clks, A] = sh_expand(R, L_max,t,p)
%% t, p are column vectors of equal length
%% The expansion is most accurate in case the pairs t and p are distributed
%% in a uniform configuration on the sphere.
global S V U invS
%%%%%%%%%%%%%%%%%%%%%% Solve Linear Least squares using SVD
warning off MATLAB:divideByZero;
[L, K ] = indices_gen(1:(L_max + 1)^2); 
M = length(L); %% number of functions in expansion
N = length(R); %% number of data points
if nargin<5, recalc = 1;end
if recalc,
    A  = zeros(N, M, 'single');
    for S = 1:length(L)
        A(:, S) = ylk_cos_sin_old(L(S),K(S),p,t)';
    end
    [U, S, V] = svd(A, 0);
    invS = 1./(S);invS(invS==inf) = 0;
end

%%% use the svd
clks = (V*invS) * (U'*R');

%%% use QR decomposition
% R = single(R);
% A = single(A);
% clks = A\R';

%%% use LU factorization
% R = single(R');
% A = single(A);
% clks = U\(L\(R(lup,:)));













