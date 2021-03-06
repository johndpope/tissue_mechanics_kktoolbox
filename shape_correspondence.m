function r = shape_correspondence( X1, X2 )
% Calculate shape correspondence like a nonlinear regression coefficient

%% make sure that both have the same length
X1 = X1(:);X2 = X2(:);
if length(X1)>length(X2),
    L_max = sqrt(length(X1)/3)-1;
    X2 = tr(X2,L_max);  %% pads with zeros
else
    L_max = sqrt(length(X2)/3)-1;
    X1 = tr(X1,L_max);  %% pads with zeros
end
%% make trs invariant
disp('Calculating transformation invariant shape 1');X1 = trs_invariance_gen(X1);
disp('Calculating transformation invariant shape 2');X2 = trs_invariance_gen(X2);
%% calculate the regression coefficient value
r = 1- sum((X1-X2).^2)/sum(X1.^2);
