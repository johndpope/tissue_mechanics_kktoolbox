classdef surface_mesh_collection
    properties
        M;
        %%% display parameters
        use_camorbit = 1;
    end
    methods
        function obj = surface_mesh_collection
           M = {};
        end
        function obj = add(obj, X, F)
            if nargin == 2,
                obj.M{end+1} = X;
            elseif nargin ==3
            obj.M{end+1} = surface_mesh(X,F);
            end
        end
        function obj = import_all_X(obj,all_X, all_F)
            for ix = 1:length(all_X)
                obj = add(obj,surface_mesh(all_X{ix}, all_F{ix}));
            end
        end
        function plot(obj)
            figure;
            rand('seed',99);
            for ix = 1:length(obj.M)
                plot_fast(obj.M{ix}, rand(1,3));hold on;
            end
            hold off;axis equal; lighting gouraud; view(3);axis off; axis vis3d;camlight;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
        end
    end
end