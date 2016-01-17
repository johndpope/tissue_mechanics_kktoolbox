function Y_LK = get_basis(t,p,gdim,L_max)
% This is a low level function used to check whether a basis exists on the hard disk to be loaded or 
% whether we have to calculate a new one (and save it).
%precalc_path = '/home/khairy/KhaledData/mwork/sh_precalc';  % update this on every system once
precalc_path = '.';

verbose = 0;
[L, K] = indices_gen(1:(L_max + 1)^2); M = length(L);N = length(t);
%str = sprintf('NOCS_gdim_%d_lmax_%d.mat', gdim, L_max);
if L_max <0 || L_max > 40, warning('L_max may be out of bounds');disp('L_max');disp(L_max);end
str = sprintf('%s/BOSH_gdim_%d_lmax_%d.mat', precalc_path, length(t), L_max);
if exist(str)==2, load(str);if verbose, disp('Loading precalculated SH basis functions');end;end;
if ~exist('Y_LK'),
    if verbose,disp('Calculating basis set functions for the first time.');end
    Y_LK  = zeros(N, M, 'single');
    if verbose, h = waitbar(0,'Setting up basis set...');end
    for S = 1:length(L),
%         if ~mod(S,50),disp([num2str(S) ' of ' num2str(length(L))]);end;
        Y_LK(:,S) = ylk_bosh(L(S),K(S),p',t')'; % usese nocs version
%        Y_LK(:,S) = ylk_cos_sin_old(L(S),K(S),p',t')'; % usese nocs version
        if verbose, waitbar(S/length(L), h);end
    end;
    if verbose, close(h);end
    save(str,'Y_LK');
elseif any(size(Y_LK)~=[N M]),
    if verbose, disp('Recalculating basis set functions to accomodate new truncation/dimension.');end
    if size(Y_LK,2)>M && size(Y_LK,1)==N, 
        disp('Truncating precalculated (large) set.');
        Y_LK = Y_LK(:,1:M);
    else,     
        Y_LK  = zeros(N, M, 'single');
        for S = 1:length(L), 
            Y_LK(:,S) = ylk_cos_sin(L(S),K(S),p',t')';
            save(str,'Y_LK');
        end
    end
else
    if verbose,disp('No need to recalculate SH basis. Using precalculated set');end
end;