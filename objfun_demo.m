function f = objfun_demo(parameters,addarg)

% OBJFUN_DEMO is a demo objective function.
%    F = OBJFUN_DEMO(P,ADDARG) returns the value of the objective
%    function. P should be a 1 x 4 or 4 x 1 vector. ADDARG is a cell
%    array that gives additional arguments passed to the function.
%    Refer to the example how ADDARG is specified.
%
%    example:
%       trueval = [2 4 6 8];
%       coef = [1 2 3 4];
%       addarg = {trueval,coef};
%       f = objfun_demo([0.5,3,2.2,5],addarg)
%

p1 = parameters(1);
p2 = parameters(2);
p3 = parameters(3);
p4 = parameters(4);

foo = addarg{1};
minp1 = foo(1);
minp2 = foo(2);
minp3 = foo(3);
minp4 = foo(4);

foo2 = addarg{2};
coef1 = foo2(1);
coef2 = foo2(2);
coef3 = foo2(3);
coef4 = foo2(4);

f = coef1.*(p1-minp1).^2 + coef2.*(p2-minp2).^2 + ...
    coef3.*(p3-minp3).^2 + coef4.*(p4-minp4).^2;
