function Xmesh = shape_mesh(xclks, yclks, zclks, gdim, plotflag)
%% Creates the mesh of the shape that will be used to 
%% generate the synthetic volume data.

L_max = round(sqrt(length(xclks))-1);
%%% preparation of quantities needed for the shape synthesis
% global Y_LK P_LK Y_LK_phi Y_LK_theta P_LK_T Y_PP Y_TT Y_TP wp wt
% gdim = 100;
% [X] = partsphere(gdim);
% whos X
% x = X(1,:);y = X(2,:); z = X(3,:); [t, p, r]= kk_cart2sph(x,y,z);

[t wt]                  = gaussquad(gdim, 0, pi);
[p wp]                  = gaussquad(gdim,0,2*pi);

[p t]                   = meshgrid(p,t);% [wp wt]                 = meshgrid(wp, wt);
[Y_LK P_LK]				= precalc_ylk_cos_sin(p, t, L_max);

c = xclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
x = sum(c.*Y_LK,3);%xp =(sum(c.*Y_LK_phi,3));xt = (sum(c.*Y_LK_theta,3));xpp = (sum(c.*Y_PP,3));xtt = (sum(c.*Y_TT,3));xtp = (sum(c.*Y_TP,3));
c = yclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
y = sum(c.*Y_LK,3);%yp =(sum(c.*Y_LK_phi,3));yt = (sum(c.*Y_LK_theta,3));ypp = (sum(c.*Y_PP,3));ytt = (sum(c.*Y_TT,3));ytp = (sum(c.*Y_TP,3));
c = zclks';h = length(c);c = c(ones(gdim^2,1),:);c = reshape(c,gdim, gdim, h); % 
z = sum(c.*Y_LK,3);%zp =(sum(c.*Y_LK_phi,3));zt = (sum(c.*Y_LK_theta,3));zpp = (sum(c.*Y_PP,3));ztt = (sum(c.*Y_TT,3));ztp = (sum(c.*Y_TP,3));
Xmesh = [reshape(x,gdim^2,1) reshape(y,gdim^2,1) reshape(z,gdim^2,1)];
if plotflag
    if     plotflag ==1 surf(x,y,z,sqrt(x.^2+y.^2+z.^2));grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    elseif plotflag ==2, surf(x,y,z,t);grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    elseif plotflag ==3, surf(x,y,z,p);grid off;daspect([1 1 1]);hold off;rotate3d; drawnow;
    end
end