function [Xsat, A, V] = satisfy_constraints_vec(clks,A_o,Y_LKtri,v_o, lambda, max_try, Areatol, Voltol, miu)
%%%% Calculate the Areas and Volumes for a large number of shapes in vectorized fashion
%%%% to assess whether they satisfy the constraints of area and volume or not.
global uC1 uC2 uC3 vC1 vC2 vC3 wC1 wC2 wC3 vec1 vec2 v1x v1y v1z v2x v2y v2z v3x v3y v3z crossqpr_vec twoA_vec
global Y_LK_C1 Y_LK_C2 Y_LK_C3 Y_LK_C1_x Y_LK_C2_x Y_LK_C3_x Y_LK_C1_y Y_LK_C2_y Y_LK_C3_y
% ok = 0;Areatol = 1;Voltol  = 0.05;
Xsat = [];
clks = clks(:)';
Xtry = zeros(max_try,length(clks));

nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
exp_vec = exp(-[1:nc]./miu);exp_vec = [exp_vec exp_vec exp_vec]/max(exp_vec);
exp_mx = exp_vec(ones(max_try,1),:);
% plot(exp_vec);drawnow;

rfact_mx = (randn(size(Xtry))-0.5).*exp_mx;

% rfact_mx = (rand(size(Xtry))-0.5);%.*(clks(ones(max_try,1),:)/max(abs(clks))) ;
Xtry = clks(ones(max_try,1),:)+lambda.*rfact_mx;   % generate large number of random shapes

nc = round(length(clks)/3);xclks = Xtry(:,1:nc)';yclks = Xtry(:,nc+1:2*nc)'; zclks = Xtry(:,2*nc+1:3*nc)';

uC1(:,1,:) = Y_LK_C1_x*xclks;
uC2(:,1,:) = Y_LK_C2_x*xclks;
uC3(:,1,:) = Y_LK_C3_x*xclks;

vC1(:,1,:) = Y_LK_C1_y*yclks;
vC2(:,1,:) = Y_LK_C2_y*yclks;
vC3(:,1,:) = Y_LK_C3_y*yclks;

wC1(:,1,:) = Y_LK_C1*zclks;
wC2(:,1,:) = Y_LK_C2*zclks;
wC3(:,1,:) = Y_LK_C3*zclks;

v1x(:,1,:) = uC2-uC1;
v1y(:,1,:) = vC2-vC1;
v1z(:,1,:) = wC2-wC1;

v2x(:,1,:) = uC3-uC1;
v2y(:,1,:) = vC3-vC1;
v2z(:,1,:) = wC3-wC1;

crossqpr_vec(:,:,:) = [(v1y.*v2z - v2y.*v1z) (v1z.*v2x - v1x.*v2z) (v1x.*v2y - v2x.*v1y)];

twoA_vec(:,1,:) = (sqrt(sum(crossqpr_vec.*crossqpr_vec,2)));   %
A = sum(squeeze(twoA_vec),1)/2';A = A'; % this is the total area
n = crossqpr_vec./twoA_vec(:,ones(3,1),:);
V = -squeeze(sum(1/3*(dot(n,[uC1 vC1 wC1], 2).*twoA_vec./2)));

Vo = 4/3*pi*(A(:)/4/pi).^(3/2);
v_red_o = v_o./(4/3*pi*(A_o/4/pi).^(3/2));
% Xsat = Xtry(abs(abs(v_red)-v_o)<Voltol&abs(abs(A)-A_o)<Areatol,:);
Xsat = Xtry(abs(abs(V)./Vo-v_red_o)<Voltol, :);% & abs(abs(A)./A_o-1)<Areatol,:);

A = abs(A)./A_o;
V = abs(V)./Vo;





