function [P] = qperimeter(X,F)
% returns the perimeters of the quad elements
if size(F,2)~=4,error('This function requires a quad mesh');end
P = 0;
P = P + sqrt(sum((X(F(:,1),:)-X(F(:,2),:)).^2,2));
P = P + sqrt(sum((X(F(:,2),:)-X(F(:,3),:)).^2,2));
P = P + sqrt(sum((X(F(:,3),:)-X(F(:,4),:)).^2,2));
P = P + sqrt(sum((X(F(:,4),:)-X(F(:,1),:)).^2,2));