classdef tfp_basis < handle
    properties
        basis = 'cos sin';
        M = 5;
        N = 5;
        gdim = 60;
        u;
        v;
        w;
        % the basis functions
        cosnu;
        sinnu;
        cosmv;
        sinmv;
        % first derivatives of basis functions
        ducosnu;
        dusinnu;
        dvcosmv;
        dvsinmv;
        % second derivatives of basis funcions
        duducosnu;
        dudusinnu;
        dvdvcosmv;
        dvdvsinmv;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Static)
    end
    
    methods
        function obj = tfp_basis(M,N, gdim)
            obj.M = M;
            obj.N = N;
            obj.gdim  = gdim;
            [obj.u wu]                      = sh_basis.gaussquad(gdim, 0, 2*pi );
            [obj.v wv]                      = sh_basis.gaussquad(gdim,0, 2*pi);
            [obj.v obj.u]                   = meshgrid(obj.v,obj.u);
            [wv wu]                         = meshgrid(wv, wu);
            obj.w                           = wv(:).*wu(:);
            
            %%% generate basis
            obj.cosnu = zeros(length(obj.u(:)), N);
            obj.sinnu = zeros(length(obj.u(:)), N);
            obj.cosmv = zeros(length(obj.v(:)), M);
            obj.sinmv = zeros(length(obj.v(:)), M);
            
            for nix = 1:N, n = nix-1;
                obj.cosnu(:,nix) = cos(n*obj.u(:));
                obj.sinnu(:,nix) = sin(n*obj.u(:));
                
                obj.ducosnu(:,nix) = -n*sin(n*obj.u(:));
                obj.dusinnu(:,nix) = n*cos(n*obj.u(:));
                
                obj.duducosnu(:,nix) = -n*n*cos(n*obj.u(:));
                obj.dudusinnu(:,nix) = -n*n*sin(n*obj.u(:));
            end
            for mix = 1:M, m = mix-1;
                obj.cosmv(:,mix) = cos(m*obj.v(:));
                obj.sinmv(:,mix) = sin(m*obj.v(:));
                
                obj.dvcosmv(:,mix) = -m*sin(m*obj.v(:));
                obj.dvsinmv(:,mix) = m*cos(m*obj.v(:));
                
                obj.dvdvcosmv(:,mix) = -m*m*cos(m*obj.v(:));
                obj.dvdvsinmv(:,mix) = -m*m*sin(m*obj.v(:));
            end
        end
    end
end