classdef shell_mesh
    properties
        m = surface_mesh; % midplane surface mesh
        u = surface_mesh; % undeformed surface mesh
        F = [];
        shell_type = 'Kirchhoff';
        Young = 1000;   % in Pa
        Poiss = 0.3; 
        uX_o = [];  % undeformed shape
    end
    methods
        function obj = shell_mesh(srf_m, srf_u, srf_g)
           % instantiate a shell_mesh object
           % srf_m ---> surface
        end
            
        
    end
end