function S_tr = nocs2bosh(S)
%%% convert the "new" CLKs to the old ones
%%% nocs2cs stands for conversion from coefficients that were calculated
%%% whithout the use of the Condon-Shortly phase factor to coefficients
%%% with. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(size(S,2),3) ~= 0,
    S = S(:,1:end-1);
end

if size(S,2) ==1,
    S_tr = nocs2cs_single(S);
else
    S_tr = zeros(size(S));
    for ix = 1:size(S,1),       % loop over the shapes
        S_tr(ix,:) = nocs2cs_single(S(ix,:));
    end
end


%%
function Xnew = nocs2cs_single(X_o)


[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);

counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        NLK_old = N_LK(L,K);
        NLK = N_LK_bosh(L,K);
        fac = NLK_old/NLK;
        fac = 1/fac;
        xclks(counter) = xclks(counter) * 1/fac;
        yclks(counter) = yclks(counter) * 1/fac;
        zclks(counter) = zclks(counter) * 1/fac;
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];