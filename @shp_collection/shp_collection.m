classdef shp_collection
    properties
        T
        use_camorbit = 1;
    end
    methods
        function obj = shp_collection
            T = {};
        end
        function obj = add(obj, s)
            obj.T{end+1} = s;
        end
        function obj = import_T(obj,T)
            L_max = shp_surface.get_L_max(T(1,:));
            basis = sh_basis(L_max, 60);
            for ix = 1:size(T,1)
                X_o = T(ix,:);
                s = shp_surface(L_max, basis);
                s.X_o = X_o;
                obj = add(obj,s);
            end
        end
        function obj = import_T_nocs(obj,T)
            L_max = shp_surface.get_L_max(T(1,:));
            basis = sh_basis(L_max, 60);
            for ix = 1:size(T,1)
                X_o = T(ix,:);
                s = shp_surface(L_max, basis);
                s.X_o = shp_surface.nocs2bosh(X_o);
                obj = add(obj,s);
            end
        end
        function obj = import_T_cs(obj,T)
            L_max = shp_surface.get_L_max(T(1,:));
            basis = sh_basis(L_max, 60);
            for ix = 1:size(T,1)
                X_o = T(ix,:);
                s = shp_surface(L_max, basis);
                s.X_o = shp_surface.nocs2bosh(shp_surface.cs2nocs(X_o));
                obj = add(obj,s);
            end
        end
        function plot(obj)
            figure;
            for ix = 1:length(obj.T)
                obj.T{ix}.use_camorbit = 0;
                plot(obj.T{ix});hold on;
            end
            hold off;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
        end
    end
end