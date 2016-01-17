function [A n d] = tri_prop(V1, V2, V3)
%%% calculate the triangle properties for a single triangle
%%% output: A: area
%%%         n: surface normal
%%%         d: total perimeter length (if required)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if nargin >1, X =[V1(:)';V2(:)';V3(:)'];
% % else
% %     X = V1;
% % end

x1 = V1(1); y1 = V1(2);z1 = V1(3);
x2 = V2(1); y2 = V2(2);z2 = V2(3);
x3 = V3(1); y3 = V3(2);z3 = V3(3);

q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2);
twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
A = sum(twoA)/2;                    % this is the total area

n = crossqpr./twoA(:,ones(3,1));
%%%
if nargout>2,           % then we also need to calculat the length of the perimeter
   d = 0;
   d1 = sqrt(sum(  (V1(:)-V2(:)).^2  ));
   d2 = sqrt(sum(  (V2(:)-V3(:)).^2  ));
   d3 = sqrt(sum(  (V1(:)-V3(:)).^2  ));
   d = [d1 d2 d3];
end


