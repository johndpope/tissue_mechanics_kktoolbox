function NMSE21 = mean_square_error_normalized(X1,X2)
% returns the normalized mean square error
MSE21 = mean_square_error_absolute(X1, X2);
MSE10 = mean(abs(X1(:)-mean(X1(:))).^2);
MSE20 = mean(abs(X2(:)-mean(X2(:))).^2);

NMSE21 = MSE21/sqrt(MSE10*MSE20);