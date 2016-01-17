%% fuse one pair of images
I1 = read_tif_stack('Y:\Shared\Processing\Ultrascope 1\11-01-15\CHN00.tif');
I2 = read_tif_stack('Y:\Shared\Processing\Ultrascope 1\11-01-15\CHN01.tif');
I_fused = image_3D_fuse(I1, I2);
write_tif_stack(uint16(mat2gray(I_fused)*65535),'Y:\Shared\Processing\Ultrascope 1\11-01-15\fused.tif');