function [RAWp] = pad_RAW(RAW,target)
%let's pad RAW


xdiff = target-size(RAW,1);ydiff = target-size(RAW,2);zdiff = target-size(RAW,4);
if mod(xdiff,2)==0, RAW = padarray(RAW,[xdiff/2 0 0 0]);else RAW = padarray(RAW, [(xdiff-1)/2 0 0 0]);RAW = padarray(RAW, [1 0 0 0], 0, 'post');end
if mod(ydiff,2)==0, RAW = padarray(RAW,[0 ydiff/2 0 0]);else RAW = padarray(RAW, [0 (ydiff-1)/2 0 0]);RAW = padarray(RAW, [0 1 0 0], 0, 'post');end
if mod(zdiff,2)==0, RAW = padarray(RAW,[0 0 0 zdiff/2]);else RAW = padarray(RAW, [0 0 0 (zdiff-1)/2]);RAW = padarray(RAW, [0 0 0 1], 0, 'post');end

RAWp = RAW;
disp(size(RAW));