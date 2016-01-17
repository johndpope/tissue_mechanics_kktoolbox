function [X,C] = bijective_newton_steps( X,C,maxiter)
% Starts with the surface mesh and uses steepest descent steps to minimize
% the local energy at vertices. The energy term may be a combination of
% bending energy, area*area dilation modulus, shear modulus*shear and other
% constraints.
% Minimizes the energy iteratively
% calculating new positions x, y,and z for the vertices.
% The algorithm uses Newton's method iteratively
clc;
seps = .05;%sqrt(eps);
gamma = 0;      % balance Eb+gamma*dAi
%% Precalculate the Jacobian sparsity pattern
h = waitbar(0,'Setting up Jacobian sparsity pattern...');
if nargin<6,filename = 'untitled';end
nvert = length(X);nfaces = length(C);
JacPat = sparse(triangulation2adjacency(C));
JacPat = sparse([JacPat JacPat JacPat]);

[i,j] = find(JacPat);
indJ = sub2ind(size(JacPat),i,j);   %indices of elements that will change
JacInd = JacPat;
JacInd(indJ) = indJ;    % matrix filled with the indices of the nonzero elements in the pattern
%indVertVal = zeros(length(indJ),3);
%pos_vec = zeros(size(indJ));        % gives the position of the variable relative to the face definition

for ix = 1:length(indJ),        % loop over the Jacobian Pattern non-zeros
    if ~mod(ix,1000),disp([num2str(ix) ' of ' num2str(length(indJ))]);end;
    [row,col] = ind2sub(size(JacPat),indJ(ix));
    %indvec = JacInd(row,:);indvec = indvec(find(indvec));
    %pos_vec(ix) = find(indvec==indJ(ix));
    %indVertVal(ix,:) = JacInd(indvec);
    indcol(ix)   = col;       %%
    indrow(ix)   = row;
    waitbar(ix/length(indJ),h);
end
close(h);
% Cv = zeros(size(indVertVal));   % the actual coordinates of the faces are stored in here
% indCv = sub2ind(size(Cv),1:length(indJ),pos_vec');
J = sparse(size(JacPat));


%% Prepare the quantities needed for the curvature calcuation
prep_curvature_calc;
X = [X(:,1);X(:,2);X(:,3)];
%% %%%%%%%%%%%%%%% Begin iterations %%%%%%%%%%%%%%%%%%%%%
%  we are taking Newton steps as follows:
%  assuming that our system of equations (constraints) is E(C), we need to find 
%  dC which makes E(C+dC) = 0;
%  According to the first order Taylor expansion:
%  Areas(C+dC) =  Areas(C) + JT(C).dC      where JT is the transpose of the Jacobian J
%  We can satisfy the system of equations by taking Newton steps:
%  We calculate dC by solving the linear system: JT(C).dC = -E(C), and then
%  incrementing C by dC.
%  The problem is that in general there are more unknowns than there are contraints (equations)
%  this means that we have many solutions. We restrict dC to the column space of E
%  to take the shortes dC possible. So dC can be written as dC = J.dS, and dS is obtianed
%  by solving the system JT(C).dC = JT(C).J(C).dS = -E(C), 
%  then we increment C by alpha J.dS
%  Note the symbol C has been subsituted by X in the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
for iter = 1:maxiter,
    disp(iter);
    
%%   %% % we need the vertex curvatures and areas at the current configuration X
    u(:) = [ X(1:(nvert))];v(:) = [X(nvert+1:2*nvert)];w(:) = [X(2*nvert+1:end)];
    crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
    twoA = sqrt(sum(crossqpr.*crossqpr,2)); A = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
    V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;
    % %%% Calculate H, dA and M (as in F.J. thesis): This requires the the sum of lengths of edges (j) for a particular vertex (i); Lij
    HidAi_mx(LVE_ix) = sqrt((u(VE_ix)-u(V_e_ix)).^2 + (v(VE_ix)-v(V_e_ix)).^2+ (w(VE_ix)-w(V_e_ix)).^2).*...   %Lij
        real(acos(n(tr1_ix,1).*n(tr2_ix,1) + n(tr1_ix,2).*n(tr2_ix,2)+ n(tr1_ix,3).*n(tr2_ix,3))).*...    %thetaij
        (sign(n(tr1_ix,1).*u(V_far_ix) + n(tr1_ix,2).*v(V_far_ix) + n(tr1_ix,3).*w(V_far_ix) - ...
        n(tr1_ix,1).*u(VE_ix) - n(tr1_ix,2).*v(VE_ix) - n(tr1_ix,3).*w(VE_ix)));
    HidAi(:) = 1/4 .* sum(HidAi_mx,2);dAij(LVT_ix) = F_areas(VT_ix)/3;
    dAi(:) = sum(dAij,2);
    Hi = HidAi./dAi;
    Energy = Hi;% + gamma*dAi;      % energies for all individual vertices	
    
    JE = JacPat;% JE(indJ) = E(indcol);  % fill in 
    for ix = 1:length(indJ),JE(indJ(ix)) = Energy(indrow(ix));end
    %%% we need the Jacobian (J) at the current step. 
%     VAL = JacPat;
%     VAL(indJ) = X(indcol);         %% Values of the current configuration
%     Cv(:) =full(VAL(indVertVal(:)));            % fill in the values without change first
%     CHG_vec = seps*X;                           %% fill in the finite difference values
%     Cv(indCv) = X(indcol)+CHG_vec(indcol);      %% now let's change the configuration for the forward differences
    CHG_vec = seps*X;
    Cv = X+seps*X;
    %% %% With the new values of u v w, Calculate H, dA and M (as in F.J. thesis): This requires the the sum of lengths of edges (j) for a particular vertex (i); Lij
    u(:) = [ Cv(1:(nvert))];v(:) = [Cv(nvert+1:2*nvert)];w(:) = [Cv(2*nvert+1:end)];
    crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
    twoA = sqrt(sum(crossqpr.*crossqpr,2)); A = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
    V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;
    HidAi_mx(LVE_ix) = sqrt((u(VE_ix)-u(V_e_ix)).^2 + (v(VE_ix)-v(V_e_ix)).^2+ (w(VE_ix)-w(V_e_ix)).^2).*...   %Lij
        real(acos(n(tr1_ix,1).*n(tr2_ix,1) + n(tr1_ix,2).*n(tr2_ix,2)+ n(tr1_ix,3).*n(tr2_ix,3))).*...    %thetaij
        (sign(n(tr1_ix,1).*u(V_far_ix) + n(tr1_ix,2).*v(V_far_ix) + n(tr1_ix,3).*w(V_far_ix) - ...
        n(tr1_ix,1).*u(VE_ix) - n(tr1_ix,2).*v(VE_ix) - n(tr1_ix,3).*w(VE_ix)));
    HidAi(:) = 1/4 .* sum(HidAi_mx,2);dAij(LVT_ix) = F_areas(VT_ix)/3;
    dAi(:) = sum(dAij,2);
    Hi_plus = HidAi./dAi;
    Energy_plus = Hi_plus;% + gamma*dAi;      % energy at vertex
    
	JE_plus = JacPat; %JE_plus(indJ) = E_plus(indcol);  % fill in 
    for ix = 1:length(indJ),JE_plus(indJ(ix)) = Energy_plus(indrow(ix));end
%     CHG = JacPat;CHG(indJ) = CHG_vec(indcol);
    J = (JE_plus-JE);
    for ix = 1:length(indJ), J(indJ(ix)) = J(indJ(ix))./CHG_vec(indrow(ix));end
    %     %%%% Calculate the Areas with the new values
%     u1 = cos(pi/2-Cv(:,1)).*cos(Cv(:,4));u2 = cos(pi/2-Cv(:,2)).*cos(Cv(:,5));u3 = cos(pi/2-Cv(:,3)).*cos(Cv(:,6));
%     v1 = cos(pi/2-Cv(:,1)).*sin(Cv(:,4));v2 = cos(pi/2-Cv(:,2)).*sin(Cv(:,5));v3 = cos(pi/2-Cv(:,3)).*sin(Cv(:,6));
%     w1 = sin(pi/2-Cv(:,1));w2 = sin(pi/2-Cv(:,2));w3 = sin(pi/2-Cv(:,3));
%     crossqpr = cross([u2-u1 v2-v1 w2-w1],[u3-u1 v3-v1 w3-w1],2);
%     areas_plus = (sqrt(sum(crossqpr.*crossqpr,2)))./2;
    %     AREASPLUS = JacPat;AREASPLUS(indJ) = areas_plus;
%     EPLUS = JacPat;EPLUS(indJ) = E_plus;
%     Jvals = (EPLUS(indJ)-E(indJ))./CHG(indJ);
    %Jvals = (E_plus-E(indrow))./CHG_vec(indcol);
    %J(indJ) = Jvals;
    %J = reshape(J,size(JacPat,1),size(JacPat,2));
    J(isnan(J)) = 0;
    J(isinf(J)) = 100;
    %%% Calculate the solution
    A  = J*J';b = -Energy;
    %[dv, flag] = bicgstab(A,b);
    %dv = bicg(A,b);
%     [L2,U2] = luinc(A,1e-3);tol = 1e-6;dv = gmres(A,b,40,tol,5,L2,U2);
    % dv = gmres(A,b,1, 5e-1, 20);
%    dv = gmres(A,b);
    %dv = pcg(A,b,1e-3,50);
    %dv = symmlq(A,b);
    tol = 1e-6;[dv] = lsqr(A,b,tol,500);
    %dv = A\b;
    %%% and take a step
    alpha = 0.1;
    X = X + alpha.*J'*dv;
    %%% plot current configuration if desired
    x = [ X(1:(nvert))];y = [X(nvert+1:2*nvert)];z = [X(2*nvert+1:end)];
    clf; patch('Vertices',[x y z],'Faces',C,'FaceVertexCData', Energy,'FaceColor','interp');
    axis square;
    %view(-123,-18);
    view(2);
    zoom(1);drawnow;

end
X = [x(:) y(:) z(:)];


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


%%    Previous calculation of Jacobian pattern. May be needed for large
%%    systems
%%%   u = zeros(length(t),1,'double');v = zeros(length(t),1,'double');
%     w = zeros(length(t),1,'double');crossqpr = zeros(size(F),'double');
% if flag
%     %%%%% Prepare initial quantities
%     JacPat = sparse(length(F),nvert);
%     for ix = 1:length(F),verts = F(ix,:);for vert = 1:length(verts),JacPat(ix,verts(vert)) = 1;end;end
%     JacPat = sparse([JacPat JacPat]);   % we have two of the same pattern
%     %indJ = find(JacPat);
%     [i,j] = find(JacPat);
%     indJ = sub2ind(size(JacPat),i,j);
%     JacInd = JacPat;JacInd(indJ) = indJ;    % fill in the jacobian pattern with the indices
%     indVertVal = zeros(length(indJ),6);
%     pos_vec = zeros(size(indJ));        % gives the position of the variable relative to the face definition
%     h = waitbar(0,'Setting up Jacobian sparsity pattern...');
%     for ix = 1:length(indJ),        % loop over the Jacobian Pattern non-zeros
%         if ~mod(ix,1000),disp([num2str(ix) ' of ' num2str(length(indJ))]);end;
%         [row,col] = ind2sub(size(JacPat),indJ(ix));
%         indvec = JacInd(row,:);indvec = indvec(find(indvec));
%         pos_vec(ix) = find(indvec==indJ(ix));
%         indVertVal(ix,:) = JacInd(indvec);
%         indcol(ix)   = col;       %%
%         indrow(ix)   = row;
%         waitbar(ix/length(indJ),h);
%     end
%     close(h);
%     Cv = zeros(size(indVertVal));   % the actual coordinates of the faces are stored in here
%     indCv = sub2ind(size(Cv),1:length(indJ),pos_vec');
%     J = sparse(size(JacPat,1), size(JacPat,2));
%     seps = sqrt(eps);
%     str = sprintf('save data_temp_newton_steps_%s J seps JacPat JacInd indJ indcol indVertVal Cv indCv indvec indrow pos_vec;',filename);eval(str);
% else
%     str = sprintf('load data_temp_newton_steps_%s ;',filename);eval(str);
% end