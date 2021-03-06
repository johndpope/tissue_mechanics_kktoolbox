function [ix1 ix2 d] = closest_point(X1,X2)
%%% determines the closes two points based on Euler distance
%%% X1 and X2 are n x 3 vectors
%%%ix1 and ix2 give the pair of closes points

C0x = X1(:,1);C0x = C0x(:,ones(size(X2,1),1));
C0y = X1(:,2);C0y = C0y(:,ones(size(X2,1),1));
C0z = X1(:,3);C0z = C0z(:,ones(size(X2,1),1));

Cdrx = X2(:,1)';Cdrx = Cdrx(ones(size(X1,1),1), :);
Cdry = X2(:,2)';Cdry = Cdry(ones(size(X1,1),1), :);
Cdrz = X2(:,3)';Cdrz = Cdrz(ones(size(X1,1),1), :);

Rd = sqrt((C0x-Cdrx).^2 + (C0y-Cdry).^2 + (C0z-Cdrz).^2);% calculate distances

ind = find(Rd==min(abs(Rd(:))));
[ix1 ix2] = ind2sub(size(Rd), ind );
d = Rd(ind);
