function H = sh4d_basis_gen_vec(Nvec,Lvec,Kvec,phi,theta, r, tau)
%% generate the basis set for the N L and K basis vectors specified and
%% return a matix H whos columns correspond to these vectors
%% Input: Nvec Lvec and Kvec must be of equal length (as obtained by
%% sh4d_indices_gen)

old = 1;
H = zeros(length(phi),length(Nvec),'single');

if old

    for ix = 1:length(Nvec)

        N = Nvec(ix);L = Lvec(ix);K = Kvec(ix);
        if N == 0 ,
            H(:,ix) = ylk_cos_sin_old(L,K,phi, theta)';
        elseif N>0
            fac = 2*pi*N/tau*r;
            H(:,ix) = ylk_cos_sin_old(L,K,phi, theta)' .*cos(fac');     %%%
        else
            fac = 2*pi*N/tau*r;
            H(:,ix) = ylk_cos_sin_old(L,K,phi, theta)' .*sin(fac');     %%%
        end
    end

else        % use the nocs version for being able to calculate derivatives

    for ix = 1:length(Nvec)

        N = Nvec(ix);L = Lvec(ix);K = Kvec(ix);
        if N == 0 ,
            H(:,ix) = ylk_cos_sin_nocs(L,K,phi, theta)';
        elseif N>0
            fac = 2*pi*N/tau*r;
            H(:,ix) = ylk_cos_sin_nocs(L,K,phi, theta)' .*cos(fac');     %%%
        else
            fac = 2*pi*N/tau*r;
            H(:,ix) = ylk_cos_sin_nocs(L,K,phi, theta)' .*sin(fac');     %%%
        end
    end

end