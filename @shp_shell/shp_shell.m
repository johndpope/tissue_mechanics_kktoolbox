classdef shp_shell < shp_surface
    properties
        F = [];
        shell_type = 'Kirchhoff';
        Young = 1000;   % in Pa
        Poiss = 0.3; 
        uX_o = [];  % undeformed shape
    end
    methods
        
    end
end