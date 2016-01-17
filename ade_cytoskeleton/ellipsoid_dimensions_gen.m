function [a, c, A, V, v] = ellipsoid_dimensions_gen(area, v_red, tol)
%%%% Calculate the a and c parameters of an ellipsoid
 %%% oblate
 if nargin == 1, v_red = 0.95;tol = 1;end
a = 0.5:.01:6.5;
c = 0.5:0.01:7.5;
[a, c] = meshgrid(a,c);
e = sqrt(1-(c.^2./a.^2));
A = 2 * pi * a.^2 + pi * c.^2./e .*log((1 + e)./(1-e));   % from mathworld for oblate ellipsoid
V = 4/3 * pi .*a.^2 .*c;
Vsph = 4/3*pi.*(A/4/pi).^(3/2);
v = V./Vsph;

v((A>(area + tol))) = nan;
v((A<(area - tol))) = nan;
v((v>(v_red+0.0001))) = nan;
v((v<(v_red-0.0001))) = nan;


indx = find(~isnan(v));
a = a(indx);c = c(indx); v = v(indx);
a = mean(a);c = mean(c);

e = sqrt(1-(c.^2./a.^2));
A = 2 * pi * a.^2 + pi * c.^2./e .*log((1 + e)./(1-e));   % from mathworld for oblate ellipsoid
V = 4/3 * pi .*a.^2 .*c;
Vsph = 4/3*pi.*(A/4/pi).^(3/2);
v = V./Vsph;

% plot3(a(:),c(:),v(:),'k*');
% xlabel('a');
% ylabel('c');
% zlabel('reduced volume');


%surf(a,c,abs(v),e);