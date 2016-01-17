function [Residual, shear, area_res] = objective(X,flag, F,F_areas_relative, DM, invDM,F_areas_o, nvert, Xfixed, AreaCalcEqual, beta)
%% Calculates the deviation of relative area from the original mesh
global pcount
global u v w %crossqpr %twoA F_areas

shear = 0;area_res = 0;
alpha = 1;if isempty(F_areas_o),beta  = 0;end;
t1 = Xfixed(1);p1 = Xfixed(2);
%t = [Xfixed(1,1); X(1:(nvert)); Xfixed(2,1)];p = [Xfixed(1,2); X(nvert+1:end); Xfixed(2,2)];
%t = [ t1;X(1:(nvert))];p = [p1;X(nvert+1:end)];
t = [ X(1:(nvert))];p = [X(nvert+1:end)];
p = mod(p,2*pi);

%[u, v, w] = kk_sph2cart(t,p ,1);% convert coordinates (on the sphere) to the Cartesian
u(:) = cos(pi/2-t(:)).*cos(p(:));
v(:) = cos(pi/2-t(:)).*sin(p(:));
w(:) = sin(pi/2-t(:));

% x1 = u(F(:,1)); y1 = v(F(:,1));z1 =  w(F(:,1));
% x2 = u(F(:,2)); y2 = v(F(:,2));z2 =  w(F(:,2));
% x3 = u(F(:,3)); y3 = v(F(:,3));z3 =  w(F(:,3));
% q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
% crossqpr = cross(q,r,2); 
% twoA = sqrt(sum(crossqpr.^2,2)); F_areas = twoA./2;Area = sum(F_areas);

crossqpr = cross([u(F(:,2))-u(F(:,1)) v(F(:,2))-v(F(:,1)) w(F(:,2))-w(F(:,1))],[u(F(:,3))-u(F(:,1)) v(F(:,3))-v(F(:,1)) w(F(:,3))-w(F(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2)); Area = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));


% n = crossqpr./twoA(:,ones(3,1));
% V = -sum(1/3*(dot(n,[u(F(:,1)) v(F(:,1)) v(F(:,1))], 2).*twoA./2));
% Vo = 4/3*pi*(Area/4/pi)^(3/2);
% RedVolume = V/Vo;       % reduced volume



if AreaCalcEqual, area_res = alpha * (F_areas);
else
    area_res = alpha * (F_areas/4/pi-F_areas_relative);
    area_res = alpha * (F_areas-F_areas_relative);
end
Residual = area_res;
%%%% Residual = alpha * (F_areas/4/pi-4*pi/length(F)); 


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
    E1 = double((-b + sqrt(b.*b - 4.* c))./2);E2 = double((-b - sqrt(b.*b - 4.* c))./2);%disp([E1 E2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_shear  = (E1 + E2 - double(ai))./(1 + double(ai));
    shear = beta_shear;
    Residual = sqrt(area_res.^2 + (beta*beta_shear).^2);     % sum the dot products
end
% Residual = sum(Residual.^2);

%%% Attach the area and volume of the sphere at the end of R
%Residual = [Residual(:);abs(Area/4/pi-1);abs(abs(RedVolume)-1)];
% absA = abs(Area/4/pi-1);
Residual = [Residual(:);0;0;0;0];


% % Plotting
if flag
    if mod(pcount,flag)==0,
%        cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');
        %figure(5);
        cla; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', F_areas_relative,'FaceColor','flat');
        view(84,-11);drawnow; axis equal; axis tight;
        
        %figure(6);plot(F_areas);hold on;plot(F_areas_relative);hold off;
        
    end
end
pcount = pcount + 1;


