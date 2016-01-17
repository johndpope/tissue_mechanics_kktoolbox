%%% test Power method for calculating the eigenvalues of a 3x3 matrix
%%%http://www.miislita.com/information-retrieval-tutorial/matrix-tutorial-3-eigenvalues-eigenvectors.html

A = [1.5 0 1;...
     -0.5 0.5 -0.5;...
     -0.5 0 0];
 x0 = [1 1 1]';
  x0 = A*x0;
 c_old = norm(Ax0);
  x0 = x0./c_old;
  
  
 