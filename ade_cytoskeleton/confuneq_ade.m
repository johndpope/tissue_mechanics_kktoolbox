function [c, ceq] = confuneq_ade(clks,A_o, q,Y_LK,v_o, da_o, lambda_shear, flag)

% global A v_red
global C itercount %A v_red
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
X = [Y_LK*xclks(:)  Y_LK*yclks(:) Y_LK*zclks(:)];
plotflag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculates the Area(A), the volume (V) and the local mean and total mean curvatures (H and h) of the shape.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global HidAi dAi HidAi_mx dAij u v w crossqpr
global V_e_ix LVE_ix VE_ix tr1_ix tr2_ix VT_ix LVT_ix V_far_ix LVE_ix_inv
%% calculation of the area and volume
u(:) = X(:,1); v(:) = X(:,2); w(:) = X(:,3);
crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],...
    [u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));   %
A = sum(twoA)/2; % this is the total area
F_areas = twoA/2;% this is the vector of face areas
n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
%Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;
tol = 2.5;
% Nonlinear inequality constraints
c = [];ix = 1;
% c(ix) =  A - A_o  - tol;ix = ix + 1;
% c(ix) = -A + A_o  - tol;ix = ix + 1;
% c(ix) =  V - v_o  - tol;ix = ix + 1;
% c(ix) = -V + v_o  - tol;ix = ix + 1;
% Nonlinear equality constraints
% ceq = [];ix = 1;ceq(ix) = abs(v_red - v_o); ix = ix + 1;
ceq = [];
ix = 1;
ceq(ix) = (V - v_o); ix = ix + 1;
ceq(ix) = (A-A_o);