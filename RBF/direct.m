function [A]=direct(cent)
%
% Form the (n+dim+1)*(n+dim+1) matrix corresponding to
% RBF interpolation with basic function phi and linears.
%
% phi is assumed given as a function of r. It is coded in
% the Matlab function phi.m
%
% Syntax [A]=direct(cent)
%
% Input
% cent n*dim array coordinates of centres
%
% Output A (n+dim+1)* Symmetric matrix of the
% (n+dim+1) linear system obtained if
% we solve for the radial
% basis function interpolant
% directly.
%
% Write the matrix A in the form
%
% B P
% A =
% P^t O
%
% where P is the polynomial bit.
%
[n dim]=size(cent);
A=zeros(n,n);
for i=1:n
for j=1:i
r=norm(cent(i,:)-cent(j,:));
temp=phi(r);
A(i,j)=temp;
A(j,i)=temp;
end
end
%
% Now the polynomial part
%
P=[ones(n,1) cent];
A = [ A P ;P' zeros(dim+1,dim+1)];