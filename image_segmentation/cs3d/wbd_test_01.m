
%filename = 'blob_smooth.tif';
filename = 'smooth_large.tif';
I = read_tif_stack(filename);
I = mat2gray(I);
%% call dgvf and gft only
disp(['Generating diffusion gradient vector field ']);
miu = 0.5;
dt = 0.009;
[un, vn, wn] = dgvf_calc(I, 100, miu, dt, 1, 1, 1, 0); % Note: distributed is set to zero here
%%
%[L num xs ys zs]  = gft_calc(un, vn, wn,0.0, 50);
%%disp(['Gradient flow tracking']);
niter = 100;
[xs,ys, zs] = meshgrid(1:size(un,2), 1:size(un,1), 1:size(un,3));
maxxr = max(xs(:));maxyr = max(ys(:));maxzr = max(zs(:));
siz = size(xs);
k = [1 cumprod(siz(1:end-1))];%Compute linear indices
for ix = 1:niter
    %disp(['Gradient flow tracking iteration:  ' num2str(uint16(ix)) '  of  '  num2str(uint16(niter))]);
    Vx = 1 + (ys(:)-1)*k(1) + (xs(:)-1)*k(2) + (zs(:)-1)*k(3);  %% this replaces sub2ind
    NV = sqrt((un(Vx).^2+vn(Vx).^2+wn(Vx).^2));
    NV(NV==0) = 1;      %% set to one to avoid nan appearing in the xs ys and zs matrices
    %%% sosi
    g = 11;
%     disp([ix round(un(Vx(g))./NV(g)) NV(g)...
%         Vx(g) xs(g) ys(g) zs(g)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    xs(:) = xs(:)+ round(un(Vx)./NV);
    xs(xs>=maxxr) = maxxr;
    xs(xs<=0) = 1;
    
    ys(:) = ys(:)+ round(vn(Vx)./NV);
    ys(ys>=maxyr) = maxyr;
    ys(ys<=0) = 1;
    
    zs(:) = zs(:)+ round(wn(Vx)./NV);
    zs(zs>=maxzr) = maxzr;
    zs(zs<=0) = 1;
    %%%
    %disp([ix sub2ind(size(xs), ys(g), xs(g), zs(g)) xs(g) ys(g) zs(g)]);
end













