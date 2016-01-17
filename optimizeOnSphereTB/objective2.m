function [Residual] = objective2(X,flag, F,F_areas_relative, DM, invDM,F_areas_o, nvert)
%% Calculates the deviation of relative area from the original mesh
global pcount
alpha = 1;
beta  = 0;
t = X(1:(nvert)); p = X(nvert+1:end); %p = mod(p,2*pi);
[u, v, w] = kk_sph2cart(t,p ,1);% convert coordinates (on the sphere) to the Cartesian
x1 = u(F(:,1)); y1 = v(F(:,1));z1 =  w(F(:,1));
x2 = u(F(:,2)); y2 = v(F(:,2));z2 =  w(F(:,2));
x3 = u(F(:,3)); y3 = v(F(:,3));z3 =  w(F(:,3));
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2); twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
F_areas = twoA./2;
%Residual = alpha * (F_areas/4/pi-F_areas_relative);     % 
% Residual = alpha * (F_areas/4/pi-4*pi/length(F)); 
 Residual = alpha * (F_areas); %
% %%% constraints
% nq = sqrt(sum(q.*q,2));nr = sqrt(sum(r.*r,2));
% ttr = acos(dot(q,r,2)./nq./nr)*180/pi;
% fac = 2;penalty = 10;
% if any(ttr>180-fac) || any(ttr<fac), disp('area violation');Residual = Residual + penalty;end
% %%%% end of constraints

if beta,        %%% distorsion by penalizing for shear

    %% MS shear energy
    global diagx1 diagy1 diagz1 SPI diag_vec_1 diag_vec diag_vec_m1 diag_vec_2 diag_vec_m2
    global diag_vec_ds diag_vec_1_ds diag_vec_m1_ds SPI2
    ai = (F_areas./F_areas_o -1);     % area_invariant
    %%%%% build the necessary transformations as in Computer Graphics p. 217.
    nfaces = length(F);V1 = [u(F(:,1)) v(F(:,1)) w(F(:,1))];V2 = [u(F(:,2)) v(F(:,2)) w(F(:,2))];V3 = [u(F(:,3)) v(F(:,3)) w(F(:,3))];
    %%% construct the translation matrix that brings V1 to the origin
    T  = speye(4*nfaces);
    % diagx1 = zeros(1,size(T,1));diagy1 = zeros(1,size(T,1));diagz1 = zeros(1,size(T,1));
    diagx1(:) = 0;diagy1(:) = 0;diagz1(:) = 0;
    diagx1(4:4:end) = -V1(:,1);T = T + spdiags(diagx1',3,size(T,1), size(T,2));
    diagy1(4:4:end) = -V1(:,2);T = T + spdiags(diagy1',2,size(T,1), size(T,2));
    diagz1(4:4:end) = -V1(:,3);T = T + spdiags(diagz1',1,size(T,1), size(T,2));
    %%% construct the necessary rotation matrix
    R = speye(4*nfaces); p1p2 = V2-V1;p1p3 = V3-V1;
    normp1p2 = sqrt(p1p2(:,1).^2 + p1p2(:,2).^2 + p1p2(:,3).^2);
    rz = p1p2./normp1p2(:,ones(3,1));
    crossp1p3p1p2 = cross(p1p3,p1p2, 2);normcrossp1p3p1p2 = sqrt(crossp1p3p1p2(:,1).^2 + crossp1p3p1p2(:,2).^2 + crossp1p3p1p2(:,3).^2);
    rx = crossp1p3p1p2./normcrossp1p3p1p2(:,ones(3,1));ry = cross(rz, rx, 2);
    diag_vec = ones(1,4 * nfaces);diag_vec(1:4:end) =rx(:,1);diag_vec(2:4:end) = ry(:,2);diag_vec(3:4:end) = rz(:,3);
    R = R + spdiags(diag_vec',0,4*nfaces, 4*nfaces);
    % diag_vec_1 = zeros(1,4 * nfaces);diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_2 = zeros(1,4 * nfaces);diag_vec_m2 = zeros(1,4 * nfaces);
    diag_vec_1(:) = 0;diag_vec_m1(:) = 0;diag_vec_2(:) = 0;diag_vec_m2(:) = 0;
    diag_vec_1(2:4:end) = rx(:,2);diag_vec_1(3:4:end) = ry(:,3);R = R + spdiags(diag_vec_1',1,size(R,1), size(R,2));
    diag_vec_m1(1:4:end) = ry(:,1);diag_vec_m1(2:4:end) = rz(:,2);R = R + spdiags(diag_vec_m1',-1,size(R,1), size(R,2));
    diag_vec_2(3:4:end) = rx(:,3);R = R + spdiags(diag_vec_2',2,size(R,1), size(R,2));
    diag_vec_m2(1:4:end) = rz(:,1);R = R + spdiags(diag_vec_m2',-2,size(R,1), size(R,2));
    R = R - SPI; %speye(4 * nfaces);
    V1pr = reshape([V1 ones(length(V1),1)]',length(V1)*4,1);
    V2pr = reshape([V2 ones(length(V2),1)]',length(V2)*4,1);
    V3pr = reshape([V3 ones(length(V3),1)]',length(V3)*4,1);
    V1r = R*T*V1pr; V2r = R * T *V2pr; V3r = R * T * V3pr; % transform
    V1r = [reshape(V1r,4,length(V1))]';V2r = [reshape(V2r,4,length(V2))]';V3r = [reshape(V3r,4,length(V3))]';
    %%% All vertices are rotated now, and all x-coordinates are zero.now we omit the x-coordinate and look at the triangles in 2-D only. note that all V1r are translated to zero.
    V2r = [V2r(:,2) V2r(:,3)]; V3r = [V3r(:,2) V3r(:,3)];
    %%% use the rotated vertices to build the necessary matrices for the calculation of the deformation
    DS = SPI2;%speye(2*nfaces);
    % diag_vec = zeros(1,2*nfaces);diag_vec_1 = zeros(1,2*nfaces);diag_vec_m1 = zeros(1,2*nfaces);
    diag_vec_ds(:) = 0;diag_vec_1_ds(:) = 0;diag_vec_m1_ds(:) = 0;
    diag_vec_ds(1:2:end) =V2r(:,1);diag_vec_ds(2:2:end) = V3r(:,2);DS = DS + spdiags(diag_vec_ds',0,2*nfaces, 2*nfaces);
    diag_vec_1_ds(2:2:end)=V3r(:,1);DS = DS + spdiags(diag_vec_1_ds',1,2*nfaces, 2*nfaces);
    diag_vec_m1_ds(1:2:end)=V2r(:,2);DS = DS + spdiags(diag_vec_m1_ds',-1,2*nfaces, 2*nfaces);
    DS = DS - SPI2;%speye(size(DS));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = 1/2 * (invDM' * DS' * DS * invDM - SPI2);;
    diagG = diag(G);eps11 = diagG(1:2:end);eps22 = diagG(2:2:end);diag1G = diag(G,1);eps12 = diag1G(1:2:end);
    b = -(eps11+eps22); c = eps11.*eps22 -eps12.*eps12; % a is always 1
    E1 = (-b + sqrt(b.*b - 4.* c))./2;E2 = (-b - sqrt(b.*b - 4.* c))./2;%disp([E1 E2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_shear  = (E1 + E2 - ai)./(1 + ai);

    Residual = Residual + beta * beta_shear;     % sum the dot products
end
 Residual = sum(Residual.^2);
% % Plotting
if flag
    if mod(pcount,flag)==0,
        cla reset;
        axis square;graphlims = [-1.1 1.1]; xlim(graphlims);ylim(graphlims); zlim(graphlims);
        patch('Vertices',[u v w],'Faces',F, 'FaceColor', 'white');
        if beta
            str = sprintf('areas pres.= %.4f   shear = %.4f',...
                sum(alpha * (F_areas/4/pi-F_areas_relative)),sum(beta * beta_shear));
            title(str);
%         else
%             title(num2str(Residual(1)));
        end
        view(-129.5, 22);  drawnow;
    end
end
pcount = pcount + 1;

