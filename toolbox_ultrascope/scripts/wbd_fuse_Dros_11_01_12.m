%%% processing \\Keller-FS\Microscopy\Ultrascope 1\11-01-08\Drosophila Histon-eGFP Embryo 2
%%% µm Bidirectional Time-Lapse 20 s-20110108-161027
clc;close all hidden; clear all
%% configure
proc_root = pwd;
deconv_out  = 'fused_deconv';
fuse_out = 'fused';
data_root = 'Z:\Ultrascope 1\11-01-12';
exp_tag = 'Drosophila_Embryo_Histone-eGFP_2µm_Bidirectional_Dual-Camera_TL_35s';
prefix = 'Drosophila_Embryo_Histone-eGFP_2µm_Bidirectional_Dual-Camera_Time-Lapse_35_s-20110111-164235';
ntime = 180;
nz    = 116;
SPM = 'SPM00'; % does not change for this experiment
SPC = 'SPC00';
ANG = 'ANG000'; % does not change for this experiment

CHN = {'CHN00', 'CHN01'};   %%% there are two chanels
CM = {'CM0', 'CM1'};
%% start processing
nproc = 1;
% if matlabpool('size') > 0, matlabpool close force local;end
% if matlabpool('size') <nproc, str = sprintf('matlabpool local %d;', nproc);eval(str);end

for tix = [123:128]% 138:142 154:170 180 ],   % loop over time points also 170 and 180
    %% generate TM string
    if      tix<10,TM = sprintf('TM0000%d',tix);
    elseif  tix<100,TM = sprintf('TM000%d',tix);
    elseif  tix<1000,TM = sprintf('TM00%d',tix);
    elseif  tix<10000,TM = sprintf('TM0%d',tix);
    end
    %% loop over the cameras
    for cmix = 1:numel(CM)
        %% make stacks and store them in I_all
        I_all = cell(2,1);
        for imix = 1:numel(CHN),
            fn = sprintf('%s\\%s\\%s\\%s\\%s\\%s-%s_%s_%s_%s_%s_PH0_%s.tif', data_root, exp_tag,SPM,TM,ANG, prefix,SPC, TM, ANG,CM{cmix},CHN{imix}, 'PLN0000');
            %im_info = imfinfo(fn);
            I = zeros([size(imread(fn)) nz]); %% initialize I
            for zix = 0:nz-1,
                %% generate PLN string
                if      zix<10,PLN = sprintf('PLN000%d',zix);
                elseif  zix<100,PLN = sprintf('PLN00%d',zix);
                elseif  zix<1000,PLN = sprintf('PLN0%d',zix);
                elseif  zix<10000,PLN = sprintf('PLN%d',zix);
                end
                I_path = sprintf('%s\\%s\\%s\\%s\\%s\\%s-%s_%s_%s_%s_%s_PH0_%s.tif', data_root, exp_tag,SPM,TM,ANG, prefix,SPC, TM, ANG,CM{cmix},CHN{imix}, PLN);
                I(:,:,zix+1) = imread(I_path,'tiff');
            end
            I_all{imix} = I;
        end
        %% fuse the stack (experiment specific)
        I_fused = image_3D_fuse(I_all{1}, I_all{2});
        fused_fn = sprintf('%s\\%s-%s_%s_%s_PH0_fused.tif', fuse_out, prefix,SPM,TM, CM{cmix});
        write_tif_stack(uint16(mat2gray(I_fused)*65535),fused_fn);
    end
end
if matlabpool('size') > 0, matlabpool close;end












