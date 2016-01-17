function [s m] = im2shp(L_max, fn)
%% call image2shp
cmd = ['C:\KK_share\mwork\kktoolbox_local\image2shp.exe ' fn];
system(cmd);
%% convert the result into shp class
[s m] = read_image2shp(L_max);

%% plot
plot(s);