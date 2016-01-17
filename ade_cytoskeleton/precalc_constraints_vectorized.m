%%%%%%%%%%%%%% Precalculate quantities for vectorized Area and Volume calcuation
global  C Y_LK uC1 uC2 uC3 vC1 vC2 vC3 wC1 wC2 wC3 u_mat v_mat w_mat vec1 vec2
global  v1x v1y v1z v2x v2y v2z v3x v3y v3z crossqpr_vec twoA_vec
uC1 = zeros(length(C),1,mxfun,'single');
uC2 = zeros(length(C),1,mxfun,'single');
uC3 = zeros(length(C),1,mxfun,'single');
vC1 = zeros(length(C),1,mxfun,'single');
vC2 = zeros(length(C),1,mxfun,'single');
vC3 = zeros(length(C),1,mxfun,'single');
wC1 = zeros(length(C),1,mxfun,'single');
wC2 = zeros(length(C),1,mxfun,'single');
wC3 = zeros(length(C),1,mxfun,'single');

v1x = zeros(length(C),1,mxfun,'single');
v1y = zeros(length(C),1,mxfun,'single');
v1z = zeros(length(C),1,mxfun,'single');
v2x = zeros(length(C),1,mxfun,'single');
v2y = zeros(length(C),1,mxfun,'single');
v2z = zeros(length(C),1,mxfun,'single');
v3x = zeros(length(C),1,mxfun,'single');
v3y = zeros(length(C),1,mxfun,'single');
v3z = zeros(length(C),1,mxfun,'single');
crossqpr_vec = zeros(length(C),3,mxfun,'single');
twoA_vec = zeros(length(C),1,mxfun,'single');

u_mat = zeros(size(Y_LK,1),mxfun, 'single');
v_mat = zeros(size(Y_LK,1),mxfun, 'single');
w_mat = zeros(size(Y_LK,1),mxfun, 'single');
vec1 = zeros(length(C), 3, mxfun,'single');
vec2 = zeros(length(C), 3, mxfun,'single');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%