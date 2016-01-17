% generate and read xml file for Ultrascope's push.exe configuration
% Author: Khaled Khairy (khairyk@janelia.hhmi.org) -- Keller lab
clc;
fn = 'push_config.xml';
%% define the struct
S = struct( 'data_root',[], 'output_root',[], 'dimensions',[],...
            'timepoint',[], 'specimen',[], 'angle', [], 'channel', [], ...
            'camera', [], 'detection_objective', [], 'structured_illumination', [], ...
            'z_step', [], 'registration', [], 'fusion', [], 'deconvolution', [], ...
            'si_reconstruction', [], 'segmentation', [], 'tracking', [], ...
            'image_compression', [], 'object_space_compression', []);
%% make some sample data
S.data_root = 'V:\temp_data\';
S.output_root = 'C:\perm_data\Experiment_date';
S.dimensions = [1024 1024 380];
S.time_point = 1;
S.specimen_name = 'MyFirstSpecimen';
S.angle = 180;
S.channel = 'GFP';
S.camera  = 'CMOS-c1';
S.detection_objective = '40x';
S.structured_illumination = 'no';
S.z_step = 10;
S.registration = 'n';
S.fusion = 'n';
S.deconvolution = 'n';
S.si_reconstruction = 'n';
S.segmentation = 'n';
S.tracking = 'n';
S.image_compression = 'y';
S.object_space_compression = 'n';
%%
write_ultrascope_xml(fn,S);
%% and read in the data
ss = read_ultrascope_xml(fn);disp(ss);