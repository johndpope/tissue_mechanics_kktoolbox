function r = kk_iseven(x)
%%% returns 0 is x is odd and 1 if x is even
if x == 0, r = 1;
else
    if mod(x,2)== 0, r = 1;
    else r = 0;
    end
end