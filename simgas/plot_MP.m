function plot_MP(X_o, facf, phys)
%% given the shape vector, in the same form as it goes into the energy function,
%%%% and the effect of the gene expression field, plot MP

Y = phys.B.Y;
Y_P = phys.B.Y_P;
Y_T = phys.B.Y_T;
Y_PP = phys.B.Y_PP;
Y_TT = phys.B.Y_TT;
Y_TP = phys.B.Y_TP;
p    = phys.B.p;
t    = phys.B.t;
w    = phys.B.w;
X_o = recover_parms(X_o,phys);
X_o = over_Basis(X_o);
[xclks yclks zclks] = get_xyz_clks(X_o);

gdimp = length(p);gdimt = length(t);
c = xclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
x = sum(c.*Y,3);

c = yclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
y = sum(c.*Y,3);

c = zclks(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
z = sum(c.*Y,3);

figure;
surf(x,y,z,reshape(facf,size(x,1), size(x,2)), 'EdgeColor', 'none');
% axis equal;axis on;axis tight;camlight;
% view(0,-90);lighting gouraud;light;camlight
