function [I, PSF] = threeD_vol_gen(X_o, x_pixel, y_pixel, z_pixel, flag, centers)
%%%%% generate a 3D volume for the shape described by X_o
%%%%% returns the 3D matix of intensities
clc;[Area,V,v_red,t,p,X_coord]=plot_sh(X_o,20);
% load X_o_pombe_cell_01_final.mat;
%  X_o = cs2nocs(X_o);
%% define global parameters
warning off;
global Iobj Iobj2 Sfft C_fft Y_LK PSF_fft x_pixel y_pixel z_pixel Xmx Xmy Xmz
global RAW ffton plotflag parm_ix full_parms avInt I_display bkg dim2 dim3
warning on;

if nargin ==1, %
    
    %%%%%%%%%%%%%%%% configure %%%%%%%%%%%%%
%     fac = 150/128;  % nucleus 02
%    fac = 120/128;  % cell_outline 01
%     fac = 150/128;  % cell_outline 02
flag = 0;fac = 10;centers = [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% fac = 10;%dim2/size(I_outline,1);
% x_pixel = 10/fac; y_pixel = 10/fac;z_pixel = 1.5;
x_pixel = fac; y_pixel = fac;z_pixel = fac;

end
%% make theoretical PSF
dim = [128 128 128];
PSF = gauss_3d_khairy(dim(1),dim(2),dim(3),2,2,6,0,5,0);% for celloutline 02
PSF_fft = fftn(PSF,[dim(1) dim(2) dim(3)]);
Iobj = zeros(dim(1),dim(2),dim(3));

%%
nc = length(X_o)/3;yclks = X_o(nc+1:2*nc);zclks = X_o(2*nc+1:end);xclks = X_o(1:nc);

%%%%%%%%%%%%%%%% configure %%%%%%%%%%%%%%%%%%%%

%xclks(1) = 0;yclks(1) = -2.3;zclks(1) = -.5;fac = 11.0; % nucleus pombe02
 xclks(1) = centers(1);yclks(1) = centers(2);zclks(1) = centers(3);fac = 1.0; % cell outline pombe02
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xclks = xclks*fac;yclks = yclks*fac;zclks = zclks*fac;

L_max = sqrt(length(xclks))-1;
% L_max = 22;
% xclks((L_max+1)^2) = 0;yclks((L_max+1)^2) = 0;zclks((L_max+1)^2) = 0;

full_parms = [xclks(:);yclks(:);zclks(:)];
nparms = (L_max+1)^2;
nc = length(xclks);
parm_ix = [1:nparms nc+1:nc+nparms 2*nc+1:2*nc+nparms];
%parm_ix  = 1:length(full_parms);%[3 13];
X_o = full_parms(parm_ix);


%% make surface
finedim = 200;
[Area,V,v_red,t,p,Xcoord,F,Y_LK] = plot_sh(xclks, yclks, zclks, finedim);axis equal
%%%%% prepare initial quantities
Xmx = zeros(length(Xcoord),1);Xmy = Xmx; Xmz = Xmx;
Iobj = zeros(size(PSF_fft));Iobj2 = Iobj;C_fft = zeros(size(PSF_fft));Sfft = zeros(size(PSF_fft));
I_display = zeros(size(PSF_fft,1), size(PSF_fft,2),3,size(PSF_fft,3));
% volume_fit_prepare_curvature_calculation; 
% load data_curvature_preparation;
%% make the first volume simulatioin
clks = full_parms; 
clks(parm_ix) = X_o;
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);


ffton = flag;
intensity_gen(xclks, yclks, zclks); % the actual volume generation

%% make RAW from the volume simulation by adding noise and calculating PSF
% for ix = 1:size(Iobj,3),
%     Iobj(:,:,ix) = Iobj(:,:,ix).*poissrnd(0.9,size(Iobj,1),size(Iobj,2));
% end;
%% use imfill to simulate a cytoplasmic dye
%Iobj = imfill(Iobj);
%%%

save Iobj_final Iobj;
I = Iobj;
I = mat2gray(I,[min(I(:)) max(I(:))]);
I = gray2ind(I,256);
% kk_montage(mat2gray(I));

%%% specific code for desired image
% Iobj = Iobj(:,:,37:111);
% 
% In = Iobj;
% load I_cell_outline.mat;
% Ic = I;
% 
% kk_montage(mat2gray(0.5*In+Ic));
% Iobj = mat2gray(Iobj);
% kk_montage(mat2gray(Iobj(:,:,17:91)));



% for ix = 1:size(Iobj,3),
%     Iobj(:,:,ix) = Iobj(:,:,ix).*poissrnd(.1,size(Iobj,1),size(Iobj,2));
% end;
% dim1 = 16;
% PSF = gauss_3d_khairy(dim1,dim1,dim1,1,1,3,0,5,0);%PSF  = mat2gray(PSF);
% dim2 = 64;dim3 = 64;
% I_outline = imresize(Iobj,dim2/size(I_outline,1));
% I_fft = fftn(I_outline,[dim2 dim2 dim3]);
% PSF_fft = fftn(PSF,[size(I_fft,1) size(I_fft,2) size(I_fft,3)]);
% RAW = ifftn(I_fft.*PSF_fft);
% I = [];I(:,:,1,:) = mat2gray(RAW);montage(mat2gray(I));
% % RAW = smooth3(RAW,'gaussian',5);
% %%
% ffton = 1;plotflag = 2;[Fun,fun_vec] = volume_fit_objective(X_o,1);
% 