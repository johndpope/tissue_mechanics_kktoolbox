function h = kk_plot3(X,c)

if nargin == 1, c = 'b';end
if strcmp(c,'r')
    h = plot3(X(:,1), X(:,2), X(:,3), 'r-');
elseif strcmp(c,'b')
    h = plot3(X(:,1), X(:,2), X(:,3), 'b-');
elseif strcmp(c,'b.')
    h = plot3(X(:,1), X(:,2), X(:,3), 'b.');
elseif strcmp(c,'k')
    h = plot3(X(:,1), X(:,2), X(:,3), 'k*');
elseif strcmp(c,'g')
    h = plot3(X(:,1), X(:,2), X(:,3), 'g*');
elseif strcmp(c,'m')
    h = plot3(X(:,1), X(:,2), X(:,3), 'm*');
elseif strcmp(c,'y')
    h = plot3(X(:,1), X(:,2), X(:,3), 'y*');
end
xlabel('x');ylabel('y');zlabel('z');axis equal;