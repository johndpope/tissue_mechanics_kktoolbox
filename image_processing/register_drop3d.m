function register_drop3d(target, source, registered, regconf, regout, drop3dpath)
% Registration using dropreg3d
% % example code preceding a call to register_drop3d
% % target = localfn{1};
% % source = localfn{2};
% % registered = [localfn{2} '_drop3d'];
% % regconf    = 'par02.txt';
% % regout     = [localfn{2} '_regout.txt'];
% % drop3dpath = '"C:\Program Files\CAMP\drop\dropreg3d.exe"';
% % I = tif2mhd(target);

% convert target and source into mhd format
I = tif2mhd(source);
%%% register using dropreg3d
str = [drop3dpath ' ' [source '.mhd'] ' ' [target '.mhd'] ' ' registered ' ' regconf ' >' regout];
system(str);
%%% convert resulting registered image to tif
mhd2tif([registered '.mhd']);
%%% delete all files that are not needed
delete([registered '.mhd'], [registered '.raw']);
delete([registered '_x.mhd'],[registered '_x.raw']);
delete([registered '_y.mhd'],[registered '_y.raw']);
delete([registered '_z.mhd'],[registered '_z.raw']);
