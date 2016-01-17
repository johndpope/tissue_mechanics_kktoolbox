function H = sh4d_basis_gen(N,L,K,phi,theta, r, tau)
%% generate the basis set

if N == 0 ,
    H = sh_basis.ylk_bosh(L,K,phi, theta);
elseif N>0
    fac = 2*pi*N/tau*r;
    H = sh_basis.ylk_bosh(L,K,phi, theta) .*cos(fac);     %%%
elseif N<0
    fac = 2*pi*(N)/tau*r;
    H = sh_basis.ylk_bosh(L,K,phi, theta) .*sin(fac);     %%%
end
