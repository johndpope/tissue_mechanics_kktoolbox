function [allX allF] = write_shps_2_mesh(S,gdim, fn)
%%% generate a wave front .obj file from a group of spherical harmonics
%%% shapes 
allX = [];
allF = [];
verbose = 0;
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape
L_max = round(sqrt(size(S,2)/3)-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});

Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
counter = 0;
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
    allX = [allX;X];
    allF = [allF;C+counter*length(X)];
    counter = counter+1;
end

% str = sprintf('%s.obj',fn);write_obj(str, allX, allF);disp('Exported obj format');
str = sprintf('%s.off',fn);write_off(str, allX, allF);disp('Exported off format');
% str = sprintf('%s.off',fn);write_off(str, allX, allF);disp('Exported off format');
% str = sprintf('%s.dxf',fn);mat2dxf(str,allX,allF);disp('Exported dxf format');
%% write to stl format
str1 = sprintf('%s.raw',fn);str2 = sprintf('%s.stl',fn);mesh2raw(str1,allX,allF);Raw2Stl(str1, str2);disp('Exported stl format');

 





