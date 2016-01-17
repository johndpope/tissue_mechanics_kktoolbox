function [Xin Xout F] = shell_fig(s, d,cflag, yval)
%% generate a standard figure assuming a shell (thickness d) from s
% s can be a file name or shp_surface object
fn = 'untitled';
if nargin ==1, d = 10;cflag = 1;yval = 10;end
if nargin ==2, cflag = 1; yval = 10;end
if nargin ==3, yval = 10;end
if strcmp('char', class(s)), 
    fn = s;
    disp('Generating shp object');
    s = shp_surface(s,5);       %% since we are not doing any calculations we use a small gaussian base point grid
    %s = truncate(s,20);
    disp('done');
else
    fn = s.name;
end
%s = truncate(s,16);
%% generate figures
dfig(1);clf;plot_pretty(s,4);
view(0,-90);camlight;
view(0,0);
myaa;drawnow;
str = sprintf('print -dtiff -r600 %s_view_01.tif',fn);eval(str);
dfig(1);view(-50,24);
myaa;drawnow;
str = sprintf('print -dtiff -r600 %s_view_02.tif',fn);eval(str);
%% generate mesh
dfig;
m = get_mesh(s,4);
[Xin Xout F]= plot_slice_y_outer(m, d,yval, cflag);
% view(-144.67, -6.56);
view(143.33, 19.4);
camlight;
lighting phong;
myaa;drawnow;
str = sprintf('print -dtiff -r600 %s_view_03.tif',fn);eval(str);
