%% generate the figure that shows the gene expression distributions of twist, snail and hkb
%%% where each is shown separately together with ftz as a reference
clc;
load d_o
load ge;
load X_o_ftz_xc.mat
X_o = d.X_o;
X_o = tr(X_o,5);
d.X_o = X_o;    % make the ftz one the shape outline


L = 25;
s = sh_surface(L,sh_basis(L,120));
sf = d.sf;
d.sf = {};
%% generate the scalar fields
s.xc = sh_surface.tr_xc(gtw,L);s = sh_rot(s,0, 0, pi);
d = add_sf(d,'twist',s);

s.xc = sh_surface.tr_xc(gsnl,L);s = sh_rot(s,0, 0, -pi/12);
d = add_sf(d,'snail',s);

s.xc = sh_surface.tr_xc(ghkb,L);
d = add_sf(d,'hkb',s);

% load X_o_ftz_xc.mat
% L = 5;
% s = sh_surface(L,sh_basis(L,120));
% s.xc = sh_surface.tr_xc(xc,L);s = sh_rot(s,0, 0, -pi);
% d = add_sf(d,'ftz',s);

az = 0;el = 0;
plot_field(d,5,1:3,1);view(az,el); % field number 4 with color red 
% plot_field(d,nico,1,2);view(az,el);saveas(gcf,'twist.fig');myaa;print -dtiff -r600 twist.tif;
