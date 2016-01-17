function S_tr = nocs2strict(S)
%%% convert the clks of shapes appearing correctly without the cs phase
%%% factor to one appearing correctly with the strict basis.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(size(S,2),3) ~= 0,
    S = S(:,1:end-1);
end

if size(S,2) ==1,
    S_tr = nocs2strict_single(S);
else
    S_tr = zeros(size(S));
    for ix = 1:size(S,1),       % loop over the shapes
        S_tr(ix,:) = nocs2strict_single(S(ix,:));
    end
end


%%
function Xnew = nocs2strict_single(X_o)

[xclks yclks zclks] = get_xyz_clks(X_o);
L_max = get_L_max(X_o);

counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
%         NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
%         NLK = N_LK(L,K);
%         fac = CS*NLK_old/NLK;
        fac = CS;
        if K<0,
            fac = CS*factorial(L-K)/factorial(L+K);     % also include CS when calculating PLK for negative K
        end

        xclks(counter) = xclks(counter) * fac;
        yclks(counter) = yclks(counter) * fac;
        zclks(counter) = zclks(counter) * fac;
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];