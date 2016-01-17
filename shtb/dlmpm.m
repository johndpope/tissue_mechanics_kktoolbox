function d = dlmpm(l, mp, m,b)
% calculate dlmpm according to Su and Coppens'94
if abs(b)<1e-10, d = 0;if mp==m, d = 1;end;end
k_sum = 0;
for k = max(0,m-mp):min(l-mp, l+m),
    k_sum = k_sum +...
            (-1)^k * factorial(l+m)/factorial(l+m-k)/factorial(k)...
            *factorial(l-m)/factorial(-m+mp +k)/factorial(l-mp-k)...
            *(cos(b/2))^(2*l-mp+m-2*k)...
            *(sin(b/2))^(2*k-m+mp);
end
d = sqrt(factorial(l+mp)*factorial(l-mp)/factorial(l+m)/factorial(l-m))...
    * (-1)^(mp-m) * k_sum;