classdef dros_embryo < shp_shell
    properties
        stage = 13;
        kb = 1;
    end
    methods
        function obj = dros_embryo(arg1,arg2, arg3)
            if nargin==0
                gdim = 120;
%                 load X_o_dros_embryo_01.mat;
                load X_o_dros_embryo_03.mat;
%                 X_o = tr(X_o,14);
                obj.L_max = shp_surface.get_L_max(X_o);
                obj.basis = sh_basis(obj.L_max, gdim);
                [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X_o);
                obj.X_o = X_o; %#ok<*PROP>
                obj = update(obj);
                obj.needs_updating = 1;
            elseif nargin ==1 && ischar(arg1)
                gdim = 120;
                load(arg1,'X_o');
                obj.L_max = shp_surface.get_L_max(X_o);
                obj.basis = sh_basis(obj.L_max, gdim);
                obj.X_o = X_o; %#ok<*PROP>
                obj = update(obj);
                obj.needs_updating = 1;
            elseif nargin==2 && strcmp(class(arg2), 'sh_basis')
                L_max = arg1;basis = arg2;
                obj.L_max = L_max;
                obj.basis = basis;
                load X_o_dros_embryo_03.mat;
                [xc yc zc] = shp_surface.get_xyz_clks(X_o);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                obj.X_o = X_o; %#ok<*PROP>
                obj = update(obj);
                obj.needs_updating = 1;
            elseif nargin==3 && strcmp(class(arg2), 'sh_basis')
                L_max = arg1;basis = arg2;X_o = arg3;
                obj.L_max = L_max;
                obj.basis = basis;
                [xc yc zc] = shp_surface.get_xyz_clks(X_o);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                obj.X_o = X_o; %#ok<*PROP>
                obj = update(obj);
                obj.needs_updating = 1;
                
            else
                error('initialize dros_embryo as >d = dros_embryo(L_max,basis)');
            end
        end
        function obj = add_sf(obj,name,s)
            if nargin ==3
                obj.sf{length(obj.sf)+1} = {name,s};
            end
        end
    end
end