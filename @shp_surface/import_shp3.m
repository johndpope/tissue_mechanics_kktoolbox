function obj = import_shp3(obj, fn, dim)
%% needs to be updated to the latest format including additional scalar fields
%% generate a shape s and read the clks from file generated using the c++ class shp_surface
if nargin ==2, dim = 20;end
fid = fopen(fn);
n_shapes = fscanf(fid,'n_shapes = %d\n',1);
L = fscanf(fid, 'L_max = %d\n', 1);
n_components = fscanf(fid,'n_components = %d\n',1);
nc = (L+1)*(L+1);
for(ix = 1:n_components)
    tags{ix} = fscanf(fid,'%s\t', 1);
end
    X = fscanf(fid,'%e\t', [n_components nc]);
fclose(fid);

%
%s = shp_surface(b);
X_o = X(1:3,:)';
s = shp_surface(L, dim);
s.X_o = X_o(:);
if n_components>3
    for ix = 4:n_components
        %disp(['Constructing scalar field ' num2str(ix-3)]);
        g = sh_surface(L,s.basis);
        g.xc = X(ix,:)';
        isf = {tags{ix}, g};
        s.sf{ix-3} = isf;
    end
end
obj = s;