function [Xsat, A, V] = satisfy_constraints_vec_ax_sym(clks,A_o,Y_LKtri,v_o, lambda, max_try, Areatol, Voltol, miu)
%%%% Calculate the Areas and Volumes for a large number of shapes in vectorized fashion
%%%% to assess whether they satisfy the constraints of area and volume or not.
global uC1 uC2 uC3 vC1 vC2 vC3 wC1 wC2 wC3 v1x v1y v1z v2x v2y v2z crossqpr_vec twoA_vec
global Y_LK_C1 Y_LK_C2 Y_LK_C3 Y_LK_C1_x Y_LK_C2_x Y_LK_C3_x Y_LK_C1_y Y_LK_C2_y Y_LK_C3_y
% ok = 0;Areatol = 1;Voltol  = 0.05;
Xsat = [];
clks = clks(:)';
Xtry = zeros(max_try,length(clks)-2);

%%% make sure that the high frequencies are penalized during the random
%%% shape generation process
nc = length(clks(3:end));
exp_vec = exp(-[1:nc]./miu);exp_vec = [exp_vec]/max(exp_vec);
exp_mx = exp_vec(ones(max_try,1),:);
rfact_mx = (randn(size(Xtry))-0.5).*exp_mx;
zclks = clks(3:end);
Xtry = zclks(ones(max_try,1),:)+lambda.*rfact_mx;   % generate large number of random shapes
xclk = clks(1)*(randn(1,max_try)+0.5);
yclk = clks(2)*(randn(1,max_try)+0.5);

%%% calculate shape properties area and volume
uC1(:,1,:) = Y_LK_C1_x*xclk;
uC2(:,1,:) = Y_LK_C2_x*xclk;
uC3(:,1,:) = Y_LK_C3_x*xclk;

vC1(:,1,:) = Y_LK_C1_y*yclk;
vC2(:,1,:) = Y_LK_C2_y*yclk;
vC3(:,1,:) = Y_LK_C3_y*yclk;

wC1(:,1,:) = Y_LK_C1*Xtry(:,1:end)';
wC2(:,1,:) = Y_LK_C2*Xtry(:,1:end)';
wC3(:,1,:) = Y_LK_C3*Xtry(:,1:end)';

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

%% test for hard constraint satisfaction
Xsat = [xclk(:) yclk(:) Xtry];
Xsat = Xsat(abs(abs(V)./Vo-v_red_o)<Voltol, :);% & abs(abs(A)./A_o-1)<Areatol,:);
A = abs(A)./A_o;
V = abs(V)./Vo;





