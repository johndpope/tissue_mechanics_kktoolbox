function MSE21 = mean_square_error_absolute(X1, X2)
MSE21 = mean(abs(X2(:)-X1(:)).^2);
