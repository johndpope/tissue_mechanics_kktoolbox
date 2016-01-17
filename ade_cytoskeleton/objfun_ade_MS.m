function [E,v_red,A, V,AreaDiff] = objfun_ade(clks, A_o, q,Y_LK,v_o, a_o_bar, lambda_shear,flag)
% Calculate the energy corresponding to the shape given by the surface harmonic expansion clks.
global C itercount %A v_red
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
X = [Y_LK*xclks(:)  Y_LK*yclks(:) Y_LK*zclks(:)];
plotflag = flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculates the Area(A), the volume (V) and the local mean and total mean curvatures (H and h) of the shape.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global HidAi dAi HidAi_mx dAij u v w crossqpr
global V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix LVE_ix_inv
%% calculation of the area and volume
u(:) = X(:,1); v(:) = X(:,2); w(:) = X(:,3);
crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2)); A = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));

Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;

% %%% Calculate H, dA and M (as in F.J. thesis): This requires the the sum of lengths of edges (j) for a particular vertex (i); Lij
HidAi_mx(LVE_ix) = sqrt((u(VE_ix)-u(V_e_ix)).^2 + (v(VE_ix)-v(V_e_ix)).^2+ (w(VE_ix)-w(V_e_ix)).^2).*...   %Lij
    real(acos(n(tr1_ix,1).*n(tr2_ix,1) + n(tr1_ix,2).*n(tr2_ix,2)+ n(tr1_ix,3).*n(tr2_ix,3))).*...    %thetaij
    (sign(n(tr1_ix,1).*u(V_far_ix) + n(tr1_ix,2).*v(V_far_ix) + n(tr1_ix,3).*w(V_far_ix) - ...
    n(tr1_ix,1).*u(VE_ix) - n(tr1_ix,2).*v(VE_ix) - n(tr1_ix,3).*w(VE_ix)));
HidAi(:) = 1/4 .* sum(HidAi_mx,2);dAij(LVT_ix) = F_areas(VT_ix)/3;dAi(:) = sum(dAij,2);
M = sum(HidAi); H = HidAi./dAi;
dA  = 2 *sum(M);             % instantaneous area difference (assuming D = 1);

% %%%% cluge  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if dA<0;
% %     disp('sign change detected');
%     dA = -dA; V = -V; v_red = -v_red; H = -H;
% end %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Calculate the energies
%%%%% constants:
kb = 2e-19;k_stretch = 5e-18/kb;k_shear = k_stretch/2;D  = 3e-3;  %units in joules and microns, then cancelled the Joules by deviding by kb
if lambda_shear,k_shear = kb/lambda_shear^2;k_shear = k_shear/kb;k_stretch = k_shear *2;end % 
%%% ADE energy
wb = 0; Eade = 0; E_shear = 0; E_stretch = 0;
Eb  = 2* sum(HidAi.^2./dAi); %bending energy as in F.J. thesis
C_o_bar = a_o_bar * q * pi/D;
EADE = q * pi/A*(dA).^2 - 2 * C_o_bar*dA;
E = Eb + EADE;% E = kb * Eb + k_bar * pi/A*(dA).^2 - 2 * C_o_bar*kb*dA;E = E/kb;%%% Full generalized bilayer couple energy (without normalizing to the area) expression as in Wortis
AreaDiff = full(D*dA);
% disp(full([Eb EADE]));

%%% MS stretch energy
%{%}
shear_fac = 1;
if shear_fac
    
global F_areas_o 
a3 = -2; a4 = 8; ai = (F_areas./F_areas_o -1);     % area_invariant
E_stretch = k_stretch/2 * sum((ai.*ai + a3 * ai .*ai.*ai +a4 * ai.*ai.*ai.*ai).*F_areas_o) ;%% multiplied by undeformed area
E = E + E_stretch*shear_fac;
% disp(E);


%% MS shear energy
global invDM DM
global diagx1 diagy1 diagz1 SPI diag_vec_1 diag_vec diag_vec_m1 diag_vec_2 diag_vec_m2
global diag_vec_ds diag_vec_1_ds diag_vec_m1_ds SPI2
%%%%% build the necessary transformations as in Computer Graphics p. 217.
nfaces = length(C);V1 = [u(C(:,1)) v(C(:,1)) w(C(:,1))];V2 = [u(C(:,2)) v(C(:,2)) w(C(:,2))];V3 = [u(C(:,3)) v(C(:,3)) w(C(:,3))];
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
beta  = (E1 + E2 - ai)./(1 + ai);b1 = 0.7; b2 = 0.75;
E_shear = k_shear *sum((beta + b1 * ai .*beta + b2 * beta.*beta).*F_areas_o) ;
E = E + shear_fac * E_shear;
% disp(E);
end
%%%%%%%%%%%%%%%%%
E = full(E);
%%%%%%%%%%%%%%%%% Plotting
if isempty(itercount),itercount = 0;end
if plotflag
    itercount = itercount + 1;
    if (mod(itercount,20)==0)||(plotflag ==2)
    cla;low = -2; high = 1; H(H<(low)) = low;    H(H>(high)) = high;
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',H, 'FaceColor', 'interp');
    str = sprintf('A:%.2f V: %.4f v: %.4f E:%.4g Eb:%.4g EADE:%.4g Estretch:%.4g Eshear:%.4g',...
    full(A), full(V),full(v_red), full(E), full(Eb),full(EADE), full(E_stretch), full(E_shear));title(str);
    axis equal;drawnow;
    end
end