clc;close all hidden; clear all
%% configure
proc_root = pwd;
fuse_out = 'D:\Khaled\Dros_11_01_12_fused_all';
data_root = 'D:\Khaled\Dros_11_01_12_fused';
exp_tag = '';
prefix = 'Drosophila_Embryo_Histone-eGFP_2µm_Bidirectional_Dual-Camera_Time-Lapse_35_s-20110111-164235';
ntime = 180;
SPM = 'SPM00'; % does not change for this experiment
SPC = 'SPC00';
ANG = 'ANG000'; % does not change for this experiment
CM = {'CM0', 'CM1'};

%%% build CM0 and CM1 pairs (names)
I_all_names = cell(ntime+1, numel(CM) + 1);
for tix = 0:ntime,   % loop over time points also 170 and 180
    %%% generate TM string
    if      tix<10,TM = sprintf('TM0000%d',tix);
    elseif  tix<100,TM = sprintf('TM000%d',tix);
    elseif  tix<1000,TM = sprintf('TM00%d',tix);
    elseif  tix<10000,TM = sprintf('TM0%d',tix);
    end
    %%% loop over the cameras
    for cmix = 1:numel(CM)
        str = sprintf('%s\\%s-%s_%s_%s_PH0_fused.tif', ...
            data_root, prefix,SPM, TM, CM{cmix});
        I_all_names{tix+1, cmix} = str;
        disp(str);
    end
     I_all_names{tix+1,end} = sprintf('%s\\%s_%s_PH0_fused.tif', fuse_out, SPM,TM);
end
%% start processing
nproc = 8;
if matlabpool('size') > 0, matlabpool close force local;end
if matlabpool('size') <nproc, str = sprintf('matlabpool local %d;', nproc);eval(str);end

%% loop over all pairs
parfor ix = 1:size(I_all_names, 1),
    %%% read the pair from disk
    I1 = read_tif_stack(I_all_names{ix,1});
    I2 = read_tif_stack(I_all_names{ix,2});
    I2 = flipdim(I2,1);
    
    %%% register the pair
    X = -50;
    Y = 286;
    tform = maketform('affine',[1 0 0; 0 1 0; X Y 1]);
    J1 = imtransform(I1,tform, 'XData',[1 size(I1,2)],'YData',[1 size(I1,1)]);
    %%% fuse the pair
    I_fused = image_3D_fuse(J1(696:1326, 244:1786,:), I2(696:1326, 244:1786,:));
    %%% write fused image to disk
    write_tif_stack(uint16(mat2gray(I_fused)*65535),I_all_names{ix,end});
end

% %% debug
% Bmx = image_max_intensity(I2, 3, 1);

%%
if matlabpool('size') > 0, matlabpool close;end












