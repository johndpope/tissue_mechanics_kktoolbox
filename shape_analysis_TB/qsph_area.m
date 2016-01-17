function sa = qsph_area(X,F)
% returns the spherical area of the quadrilateral elements in F
X = normalizeVector3d(X); % do this in vectorized fashion
a11 = zeros(size(F,1),1);
a12 = zeros(size(F,1),1);
a13 = zeros(size(F,1),1);

a21 = zeros(size(F,1),1);
a22 = zeros(size(F,1),1);
a23 = zeros(size(F,1),1);


for ix = 1:size(F,1)
    % the first triangle
    A = X(F(ix,1),:)';B = X(F(ix,2),:)';C = X(F(ix,3),:)';
    a11(ix) = sph_ang(A,B,C);
    a12(ix) = sph_ang(B,C,A);
    a13(ix) = sph_ang(C,A,B);
    
    
%     xa11 = B'*C-A'*C*A'*B; ya11 = A'*cross(B,C);% a1
%     xa12 = A'*C-B'*C*B'*A; ya12 = B'*cross(A,C);% a2
%     xa13 = B'*A-C'*A*C'*B; ya13 = C'*cross(B,A);% a3
%     a11(ix) = atan2(ya11, xa11);
%     a12(ix) = atan2(ya12, xa12);
%     a13(ix) = atan2(ya13, xa13);
    % the second triangle
    A = X(F(ix,1),:)';B = X(F(ix,4),:)';C = X(F(ix,3),:)';
    a21(ix) = sph_ang(A,B,C);
    a22(ix) = sph_ang(B,C,A);
    a23(ix) = sph_ang(C,A,B);
%     xa21 = B'*C-A'*C*A'*B; ya21 = A'*cross(B,C);% a1
%     xa22 = A'*C-B'*C*B'*A; ya22 = B'*cross(A,C);% a2
%     xa23 = B'*A-C'*A*C'*B; ya23 = C'*cross(B,A);% a3
%     a21(ix) = atan2(ya21, xa21);
%     a22(ix) = atan2(ya22, xa22);
%     a23(ix) = atan2(ya23, xa23);
    
    
    da1(ix) = a11(ix) + a12(ix) + a13(ix) -pi;
    da2(ix) = a21(ix) + a22(ix) + a23(ix) -pi;

end
sa = real(da1 + da2)';
% figure;
% plot([da1(:) da2(:) sa(:)],'-');

%%
function [alpha , a, b, c] = sph_ang(A,B,C)
% assuming A B C are already normalized

a = acos(dot(B(:)', C(:)', 2));
b = acos(dot(A(:)', C(:)', 2));
c = acos(dot(A(:)', B(:)', 2));

%%% if not normalized: use this:
%     a = vectorAngle3d(B(:)', C(:)');
%     b = vectorAngle3d(A(:)', C(:)');
%     c = vectorAngle3d(A(:)', B(:)');
    

alpha = acos((cos(a) - cos(b).*cos(c))/(sin(b).*sin(c)));









