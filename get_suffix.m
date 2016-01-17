function str = get_suffix(nz, iter)
% returns the suffix given a particular number of zeros nz, and an iterator
% iter
str = [];
if nz == 2,
    if iter<10, str = sprintf('0%d',iter);
    elseif iter<100, str = sprintf('%d',iter);
    end
end
if nz == 3,
    if iter<10, str = sprintf('00%d',iter);
    elseif iter<100, str = sprintf('0%d',iter);
    elseif iter<1000, str = sprintf('%d',iter);
    end
end
if nz == 4,
    if iter<10, str = sprintf('000%d',iter);
    elseif iter<100, str = sprintf('00%d',iter);
    elseif iter<1000, str = sprintf('0%d',iter);
    elseif iter<10000, str = sprintf('%d',iter);
    end
end