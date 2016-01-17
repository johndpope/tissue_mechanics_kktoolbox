function [alpha, beta] = tri_shear(V1, V2, V3, v1, v2, v3)
persistent invDM
%% assumes input to be vectors nx3, with each row representing a different triangle
% V1, V2, V3 are vertices of the undeformed triangles
%% returns (for every triangle)
% alpha: stretch invariant
% beta: shear invariant
nt = size(V1,1);    % number of triangles
if isempty(invDM),  % do we need to calculate invDM or not
    invDM = zeros(nt,9);
    dm = [V3-V1 V2-V1 ones(nt,3)];
    for ix = 1:nt
       invDM(ix,:) = reshape(inv(reshape(dm(ix,:),3,3)),1,9);
    end
end
DS = [v3-v1 v2-v1 ones(nt,3)];
% construct the right Cauchy-Green tensor
%C = invDM' * DS' * DS * invDM;
DSinvDM = kk_mx_mult(DS,invDM);
DSTDSinvDM = kk_mx_mult(kk_transpose(DS),DSinvDM);
C = kk_mx_mult(kk_transpose(invDM), DSTDSinvDM);
[res] = real(eig3(reshape(C',3,3,size(C,1)))); % uses vectorized  Cardan formula for root finding
res = sort(res,1)';
lambda = sqrt(res); % take square root for finding the principal stretches
%% shear and stretch invariants beta and ai
beta  = real((lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2));
alpha = real(lambda(:,1).*lambda(:,2) - 1);     % area_invariant (as function of principal stretches)

function r = kk_transpose(a)
%%% a is nx9
%%% returns nx9
r = [a(:,1) a(:,4) a(:,7) a(:,2) a(:,5) a(:,8) a(:,3) a(:,6) a(:,9)];
function r = kk_mx_mult(a,b)
%%% a and b are assumed to be nx9 arrays -- we are vectorizing matrix
%%% multiplication.
%%% returns an nx9 array
r1 = a(:,1).*b(:,1) + a(:,4).*b(:,2) + a(:,7).*b(:,3);
r4 = a(:,1).*b(:,4) + a(:,4).*b(:,5) + a(:,7).*b(:,6);
r7 = a(:,1).*b(:,7) + a(:,4).*b(:,8) + a(:,7).*b(:,9);
r2 = a(:,2).*b(:,1) + a(:,5).*b(:,2) + a(:,8).*b(:,3);
r5 = a(:,2).*b(:,4) + a(:,5).*b(:,5) + a(:,8).*b(:,6);
r8 = a(:,2).*b(:,7) + a(:,5).*b(:,8) + a(:,8).*b(:,9);
r3 = a(:,3).*b(:,1) + a(:,6).*b(:,2) + a(:,9).*b(:,3);
r6 = a(:,3).*b(:,4) + a(:,6).*b(:,5) + a(:,9).*b(:,6);
r9 = a(:,3).*b(:,7) + a(:,6).*b(:,8) + a(:,9).*b(:,9);
r = [r1 r2 r3 r4 r5 r6 r7 r8 r9];