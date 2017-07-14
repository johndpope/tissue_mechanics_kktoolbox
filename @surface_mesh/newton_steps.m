function [t,p] = newton_steps( t,p,F,stepfac, maxiter,flag,filename, verbose)
% Minimizes the areas of triangles (on the sphere) by iteratively
% calculating new positions t and p for the vertices.
% The algorithm uses Newton's method iteratively
if nargin<6,filename = 'untitled';verbose = 0;end
if nargin<7, verbose = 0;end
nvert = length(t);nfaces = length(F);p = mod(p,2*pi);X = [t(:);p(:)];
    u = zeros(length(t),1,'double');v = zeros(length(t),1,'double');
    w = zeros(length(t),1,'double');crossqpr = zeros(size(F),'double');
    if flag==0, % check whether the filename exists or not, if not set flag =1
        str = sprintf('data_temp_newton_steps_%s',filename);
        if exist(str,'file')==0, flag=1;
        end
    end
if flag
    %%%%% Prepare initial quantities
    JacPat = sparse(length(F),nvert);
    for ix = 1:length(F),verts = F(ix,:);for vert = 1:length(verts),JacPat(ix,verts(vert)) = 1;end;end
    JacPat = sparse([JacPat JacPat]);   % we have two of the same pattern
    %indJ = find(JacPat);
    [i,j] = find(JacPat);
    indJ = sub2ind(size(JacPat),i,j);
    JacInd = JacPat;JacInd(indJ) = indJ;    % fill in the jacobian pattern with the indices
    indVertVal = zeros(length(indJ),6);
    pos_vec = zeros(size(indJ));        % gives the position of the variable relative to the face definition
    if verbose, h = waitbar(0,'Setting up Jacobian sparsity pattern...');end
    indcol = zeros(length(indJ), 1, 'uint16');
    indrow = zeros(length(indJ), 1, 'uint16');
    
    %disp('Generating Jacobian sparsity pattern...');
    for ix = 1:length(indJ),        % loop over the Jacobian Pattern non-zeros
        if ~mod(ix,5000),disp([num2str(ix) ' of ' num2str(length(indJ))]);end;
        [row,col] = ind2sub(size(JacPat),indJ(ix));
        indvec = JacInd(row,:);
        indvec = indvec(indvec>0);
        pos_vec(ix) = find(indvec==indJ(ix));
        indVertVal(ix,:) = JacInd(indvec);
        indcol(ix)   = col;       %%
        indrow(ix)   = row;
        if verbose, waitbar(ix/length(indJ),h);end
    end
if verbose, close(h);end
    Cv = zeros(size(indVertVal));   % the actual coordinates of the faces are stored in here
    indCv = sub2ind(size(Cv),1:length(indJ),pos_vec');
    J = sparse(size(JacPat,1), size(JacPat,2));
    seps = sqrt(eps);
    if flag==1, str = sprintf('save data_temp_newton_steps_%s J seps JacPat JacInd indJ indcol indVertVal Cv indCv indvec indrow pos_vec;',filename);eval(str);end
else
    str = sprintf('load data_temp_newton_steps_%s ;',filename);eval(str);
    if verbose, disp('Preallocating Jakobian...');end
    J(indJ) = ones(size(indJ));     % just to allocate space beforehand
    if verbose, disp('Done !');end
end


%%%%%%%%%%%%%%%%% Begin iterations %%%%%%%%%%%%%%%%%%%%%
%  we are taking Newton steps as follows:
%  assuming that our system of equations (constraints) is Areas(S), we need to find 
%  dS which makes Areas(S+dS) = 0;
%  According to the first order Taylor expansion:
%  Areas(S+dS) =  Areas(S) + JT(S).dS      where JT is the transpose of the Jacobian J
%  We can satisfy the system of equations by taking Newton steps:
%  We calculate dS by solving the linear system: JT(S).dS = -Areas(S), and then
%  incrementing S by dS.
%  The problem is that in general there are more unknowns than there are contraints (equations)
%  this means that we have many solutions. We restrict dS to the column space of A
%  to take the shortes dS possible. So dS can be written as dS = J.dV, and dV is obtianed
%  by solving the system JT(S).dS = JT(S).J(S).dV = -Areas(S), 
%  then we increment S by J.dV
%  Note the symbol S has been subsituted by X in the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
warning off;
for iter = 1:maxiter,
    if verbose, disp(iter);end
    t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
    
    %%% we need the triangle areas (Areas) at the current configuration.
    %%% Areas is a column vector of length N-4 (where N is the # of vertices).
    u(:) = cos(pi/2-t(:)).*cos(p(:));v(:) = cos(pi/2-t(:)).*sin(p(:));w(:) = sin(pi/2-t(:));
    crossqpr = cross([u(F(:,2))-u(F(:,1)) v(F(:,2))-v(F(:,1)) w(F(:,2))-w(F(:,1))],[u(F(:,3))-u(F(:,1)) v(F(:,3))-v(F(:,1)) w(F(:,3))-w(F(:,1))],2);
    Areas = (sqrt(sum(crossqpr.*crossqpr,2)))./2;
    %%% we need the Jacobian (J) at the current step. J is a (N-4) x N matrix
    %CHG_vec = sqrt(eps)*sign(X).*(abs(X));
    CHG_vec = seps*X;                           %% fill in the finite difference values
    VAL = JacPat;VAL(indJ) = X(indcol);         %% Values of the current configuration
    Cv(:) =full(VAL(indVertVal(:)));            % fill in the values without change first
    Cv(indCv) = X(indcol)+CHG_vec(indcol);
    
    %%%% Calculate the Areas with the new values
    u1 = cos(pi/2-Cv(:,1)).*cos(Cv(:,4));u2 = cos(pi/2-Cv(:,2)).*cos(Cv(:,5));u3 = cos(pi/2-Cv(:,3)).*cos(Cv(:,6));
    v1 = cos(pi/2-Cv(:,1)).*sin(Cv(:,4));v2 = cos(pi/2-Cv(:,2)).*sin(Cv(:,5));v3 = cos(pi/2-Cv(:,3)).*sin(Cv(:,6));
    w1 = sin(pi/2-Cv(:,1));w2 = sin(pi/2-Cv(:,2));w3 = sin(pi/2-Cv(:,3));
    crossqpr = cross([u2-u1 v2-v1 w2-w1],[u3-u1 v3-v1 w3-w1],2);
    areas_plus = (sqrt(sum(crossqpr.*crossqpr,2)))./2;
    %     AREASPLUS = JacPat;AREASPLUS(indJ) = areas_plus;
    %     Jvals = (AREASPLUS(indJ)-AREAS(indJ))./CHG(indJ);
    Jvals = (areas_plus-Areas(indrow))./CHG_vec(indcol);
    J(indJ) = Jvals;
    %J = reshape(J,size(JacPat,1),size(JacPat,2));
    J(isnan(J)) = 0;
    J(isinf(J)) = 100;
    %%% Calculate the solution
    %A  = J*J';b = -Areas;
% % % %     %[dv, flag] = bicgstab(A,b);
% % % %     %dv = bicg(A,b);
% % % %     %[L2,U2] = luinc(A,1e-3);tol = 1e-6;dv = gmres(A,b,40,tol,5,L2,U2);
% % % %     %dv = gmres(A,b,1, 5e-1, 20);
        
%     [dv,flag4,relres4,iter4,resvec4]  = gmres(J*J',-Areas,2, 5e-1, 20);
    
    [dv,flag4,relres4,iter4,resvec4]  = gmres(J*J',-Areas,1, 5e0, 10);


%    dv = gmres(A,b);
    %dv = pcg(A,b,1e-3,10);
    %dv = symmlq(A,b);
    %tol = 1e-6;[dv] = lsqr(A,b,tol);
    %dv = A\b;
    %%% and take a step
    %stepfac = .5;
    X = X + stepfac.*J'*dv;

    %%% plot current configuration if desired
    t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
    u(:) = cos(pi/2-t(:)).*cos(p(:));v(:) = cos(pi/2-t(:)).*sin(p(:));w(:) = sin(pi/2-t(:));
    if verbose==2, clf; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');axis square;view(-276,4);zoom(1);drawnow;end
    

end
%%% plot current configuration if desired
t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
u(:) = cos(pi/2-t(:)).*cos(p(:));v(:) = cos(pi/2-t(:)).*sin(p(:));w(:) = sin(pi/2-t(:));
if verbose, clf; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');axis square;view(-123,-18);zoom(1);drawnow;end

warning on;

% function Jac = fdj(fun,t,p)
% %% calculate the forward difference Jacobian for a function that takes t,p as its argument
% %% and returns a vector
% f   = feval(fun,X);
% Jac = zeros(length(f),length(X));
% CHG = sqrt(eps)*sign(X).*(abs(X));
% Xplus = X;
% for ix = 1:length(X)
%     Xplus(ix) = X(ix) + CHG(ix);
%     fplus = feval(fun,Xplus);
%     Jac(:,ix) = (fplus-f)/CHG(ix);
%     Xplus(ix) = X(ix);
% end