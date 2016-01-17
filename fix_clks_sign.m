function [X1  p count indx] = fix_clks_sign(X_o, X1, thresh)
% just compare the values and if the percent difference of the absolute 
% values of corresponding coefficients lies below thresh, then flip the
% sign
indx = [];
count = 0;
if numel(X_o)~=numel(X1),error('Inputs must be same size');end
for ix = 1:numel(X_o)
    d = abs(X_o(ix))-abs(X1(ix));
    p(ix) = abs(d/X_o(ix)) * 100;
    if p(ix)<thresh
        if sign(X_o(ix))==1 && sign(X1(ix))==-1, 
            X1(ix) = -X1(ix);
            count = count + 1;
            indx = [indx;ix];
        elseif sign(X_o(ix)) == -1 && sign(X1(ix))==1, 
            X1(ix) = -X1(ix);
            count = count + 1;
            indx = [indx;ix];
        end
    end
end
disp(['Number of signs switched: ' num2str(count)]);
disp([indx X_o(indx)' X1(indx)' p(indx)']);