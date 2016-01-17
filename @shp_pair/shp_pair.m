classdef shp_pair
    properties
        S1;S2;
        n = 200;
        needs_updating = 1;
        %%% display parameters
        use_camorbit = 1;
        edge_color   = 'none';
    end
    methods(Static)
    end
    methods
        function obj = shp_pair(S1, S2)
            if nargin == 2 && strcmp('shp_surface', class(S1)) && strcmp('shp_surface', class(S2)),
                S1.use_camorbit = 0;
                S2.use_camorbit = 0;
                obj.S1 = S1;
                obj.S2 = S2;
            else
                error('shp_pair must be initialized using two shp_surface objects');
            end
        end
        function d = shape_distance(obj)
           if needs_updating,
               obj.S1trs = r_inv(obj.S1);
               obj.S2trs = r_inv(obj.S2);
               needs_updating = 0;
           end
           d = sum((obj.S1trs.X_o - obj.S2trs.X_o).^2);
        end
        function [m z t] = procrustus_distance(obj)
           if needs_updating,
               obj.S1trs = r_inv(obj.S1);
               obj.S2trs = r_inv(obj.S2);
               needs_updating = 0;
           end
            %% calculate the dissimilarity using procrustes analysis
            [m z t] = procrustes(obj.S1trs.X_o,obj.S2trs.X_o,'reflection',false);
             % Z = b*Y*T + c , i.e., z = scale*X2*rotation + translation;
        end
        function plot(obj)
            figure;
            plot(obj.S1);hold on;
            plot(obj.S2);
            hold off;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
            axis on;rotate3d;
        end
        function d = closest_distance(obj)
            %%%% function under development
            %%%% NOT DONE YET
            
            [x y z t p] = random_points_on_sphere(200);
            %%% determine L_max
            L_max = max([obj.S1.L_max obj.S2.L_max]);
            [Y, P]= sh_basis.ylk_cos_sin_bosh(p', t', L_max);
            Y = squeeze(Y);
            
            lb = (obj.S1.L_max + 1)*(obj.S1.L_max + 1);
            x1 = Y(:,1:lb)*(obj.S1.xc);
            y1 = Y(:,1:lb)*(obj.S1.yc);
            z1 = Y(:,1:lb)*(obj.S1.zc);
            
            lb = (obj.S2.L_max + 1)*(obj.S2.L_max + 1);
            x2 = Y(:,1:lb)*(obj.S2.xc);
            y2 = Y(:,1:lb)*(obj.S2.yc);
            z2 = Y(:,1:lb)*(obj.S2.zc);
            
            %%d  = min(sqrt((x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2));
        end
        function s = merge_simple(obj)
            s = obj.S1;
            X_o = (obj.S1.X_o(:) + obj.S2.X_o(:))/2;
            s.X_o = X_o;
            s = update(s);
        end
    end
end















