classdef tfp_surface
    properties
        basis;
        M;
        N;
        gdim = 60;          % default if no basis is given
        
        ax;bx;cx;dx;
        ay;by;cy;dy;
        az;bz;cz;dz;
        
        x;y;z;
        
        A;V;v;H;KG;T;h;S;Eb;Ro;m
        X_o;X_1;X_2;ang;res_p;res_o;
        xu, yu, zu, xv, yv, zv , xuu, yuu, zuu, xvv, yvv, zvv, xuv, yuv, zuv;
        needs_updating = 1;
        %%% display parameters
        use_camorbit = 1;
        edge_color   = 'none';
    end
    methods(Static)
        
    end
    methods
        function obj = tfp_surface(M, N, basis)
            if nargin == 0,         % then make a torus
                obj.M = 5;
                obj.N = 5;
                obj.gdim = 60;
                obj.basis = tfp_basis(obj.M, obj.N, obj.gdim);
            else
                obj.M = M;
                obj.N = N;
                obj.gdim = basis.gdim;
                obj.basis = basis;
            end
            % set up the coefficients that would provide a simple torus
            obj.ax = zeros(obj.N,1);obj.ay = zeros(obj.N,1);obj.az = zeros(obj.N,1);
            obj.bx = zeros(obj.N,1);obj.by = zeros(obj.N,1);obj.bz = zeros(obj.N,1);
            
            obj.cx = zeros(obj.M,1);obj.cy = zeros(obj.M,1);obj.cz = zeros(obj.M,1);
            obj.dx = zeros(obj.M,1);obj.dy = zeros(obj.M,1);obj.dz = zeros(obj.M,1);
            
            %% for x(u,v)
            fac = sqrt(2);      % sqrt(2) gives clifford torus
            A1x = fac;
            C0x = 1;
            C1x = 1/fac;
            obj.ax(2) = A1x;
            obj.cx(1) = C0x;
            obj.cx(2) = C1x;
            %% for y(u,v)
            B1y = fac;
            C0y = 1;
            C1y = 1/fac;
            obj.by(2) = B1y;
            obj.cy(1) = C0y;
            obj.cy(2) = C1y;
            %% for z(u,v)
            A0z = 1.0;
            D1z = 1;
            obj.az(1) = A0z;
            obj.dz(2) = D1z;
            
            obj.X_o = get_X_o(obj);
            obj = update(obj);
            obj.needs_updating = 0;
            
            
        end
        function obj = prepare_energy_calc(obj)
                obj.x = zeros(length(obj.basis.u(:)), 1);
                obj.y = zeros(length(obj.basis.u(:)), 1);
                obj.z = zeros(length(obj.basis.u(:)), 1);
                
                obj.xu = zeros(length(obj.basis.u(:)), 1);
                obj.yu = zeros(length(obj.basis.u(:)), 1);
                obj.zu = zeros(length(obj.basis.u(:)), 1);
                
                obj.xv = zeros(length(obj.basis.u(:)), 1);
                obj.yv = zeros(length(obj.basis.u(:)), 1);
                obj.zv = zeros(length(obj.basis.u(:)), 1);
                
                obj.xuu = zeros(length(obj.basis.u(:)), 1);
                obj.yuu = zeros(length(obj.basis.u(:)), 1);
                obj.zuu = zeros(length(obj.basis.u(:)), 1);
                
                obj.xvv = zeros(length(obj.basis.u(:)), 1);
                obj.yvv = zeros(length(obj.basis.u(:)), 1);
                obj.zvv = zeros(length(obj.basis.u(:)), 1);
                
                obj.xuv = zeros(length(obj.basis.u(:)), 1);
                obj.yuv = zeros(length(obj.basis.u(:)), 1);
                obj.zuv = zeros(length(obj.basis.u(:)), 1);
        end
        function [obj, SSn] = update_fast(obj)
                obj.x(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.y(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.z(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                obj.xu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.yu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.zu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                obj.xv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.yv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.zv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                obj.xuu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.yuu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.zuu(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                obj.xvv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.yvv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.zvv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                obj.xuv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.yuv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                obj.zuv(:) = 0;% = zeros(length(obj.basis.u(:)), 1);
                
                for mix = 1:obj.M,
                    for nix = 1:obj.N,
                        %%% calculate the coordinates
                        obj.x = obj.x + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        obj.y = obj.y + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        obj.z = obj.z + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the first derivatives w.r.t. u
                        obj.xu = obj.xu + (obj.ax(nix)*obj.basis.ducosnu(:,nix) + obj.bx(nix)*obj.basis.dusinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        obj.yu = obj.yu + (obj.ay(nix)*obj.basis.ducosnu(:,nix) + obj.by(nix)*obj.basis.dusinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        obj.zu = obj.zu + (obj.az(nix)*obj.basis.ducosnu(:,nix) + obj.bz(nix)*obj.basis.dusinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the first derivatives w.r.t. v
                        obj.xv = obj.xv + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvsinmv(:,mix));
                        obj.yv = obj.yv + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvsinmv(:,mix));
                        obj.zv = obj.zv + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvsinmv(:,mix));
                        
                        %%% calculate the second derivatives w.r.t. u
                        obj.xuu = obj.xuu + (obj.ax(nix)*obj.basis.duducosnu(:,nix) + obj.bx(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        obj.yuu = obj.yuu + (obj.ay(nix)*obj.basis.duducosnu(:,nix) + obj.by(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        obj.zuu = obj.zuu + (obj.az(nix)*obj.basis.duducosnu(:,nix) + obj.bz(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the second derivatives w.r.t. v
                        obj.xvv = obj.xvv + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvdvsinmv(:,mix));
                        obj.yvv = obj.yvv + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvdvsinmv(:,mix));
                        obj.zvv = obj.zvv + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvdvsinmv(:,mix));
                        
                        %%% calculate duv
                        obj.xuv = obj.xuv + (obj.ax(nix)*obj.basis.ducosnu(:,nix) + obj.bx(nix)*obj.basis.dusinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvsinmv(:,mix));
                        obj.yuv = obj.yuv + (obj.ay(nix)*obj.basis.ducosnu(:,nix) + obj.by(nix)*obj.basis.dusinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvsinmv(:,mix));
                        obj.zuv = obj.zuv + (obj.az(nix)*obj.basis.ducosnu(:,nix) + obj.bz(nix)*obj.basis.dusinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvsinmv(:,mix));
                        
                    end
                end
                obj.x = reshape(obj.x,size(obj.basis.u));
                obj.y = reshape(obj.y,size(obj.basis.u));
                obj.z = reshape(obj.z,size(obj.basis.u));
                
                obj.xu = reshape(obj.xu,size(obj.basis.u));
                obj.yu = reshape(obj.yu,size(obj.basis.u));
                obj.zu = reshape(obj.zu,size(obj.basis.u));
                
                obj.xv = reshape(obj.xv,size(obj.basis.u));
                obj.yv = reshape(obj.yv,size(obj.basis.u));
                obj.zv = reshape(obj.zv,size(obj.basis.u));
                
                obj.xuu = reshape(obj.xuu,size(obj.basis.u));
                obj.yuu = reshape(obj.yuu,size(obj.basis.u));
                obj.zuu = reshape(obj.zuu,size(obj.basis.u));
                
                obj.xvv = reshape(obj.xvv,size(obj.basis.u));
                obj.yvv = reshape(obj.yvv,size(obj.basis.u));
                obj.zvv = reshape(obj.zvv,size(obj.basis.u));
                
                obj.xuv = reshape(obj.xuv,size(obj.basis.u));
                obj.yuv = reshape(obj.yuv,size(obj.basis.u));
                obj.zuv = reshape(obj.zuv,size(obj.basis.u));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% Calculate first and second fundamental forms%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                X  =[obj.x(:) obj.y(:) obj.z(:)];
                Xt = [obj.xu(:) obj.yu(:) obj.zu(:)];
                Xp = [obj.xv(:) obj.yv(:) obj.zv(:)];
                Xpp = [obj.xvv(:) obj.yvv(:) obj.zvv(:)];
                Xtp = [obj.xuv(:) obj.yuv(:) obj.zuv(:)];
                Xtt = [obj.xuu(:) obj.yuu(:) obj.zuu(:)];
                
                E = dot(Xt,Xt,2);
                F = dot(Xt,Xp,2);
                G = dot(Xp,Xp,2);
                SS = (cross(Xt,Xp,2));
                SSn = sqrt(E.*G-F.*F);
                n = SS./SSn(:,ones(1,3));
                L = dot(Xtt,n,2);
                M = dot(Xtp,n,2);
                N = dot(Xpp,n,2);
                %%%%%%%%%%%%%%%%%%%%%%%% update geometrical properties
                obj.V = abs(1/3.*sum(obj.basis.w.*(dot(X,n,2)).*SSn));
                obj.A = sum(obj.basis.w.*SSn);
                H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
                KG = (L.*N - M.*M)./(E.*G-F.*F);
                %% quantities of interest
                obj.h = 1./obj.A.*sum(H(:).*obj.basis.w.*SSn);                	% normalized total mean curvature
                obj.T = sum(sum(KG(:).*obj.basis.w.*SSn))/2/pi;            		% total Gaussian curvature (constant for topology) should be zero for torus
                %obj.S = (2.*obj.H.^2-obj.KG).^(0.5);                     		% curvedness -- uncomment as needed
                %obj.Eb = 1/2*sum((2.*H(:)).^2.*obj.basis.w.*SSn)./4/pi/pi;      % bending energy relative to that of the Clifford torus
                
                r = sqrt(obj.A/4/pi);% Area of a sphere = 4 * pi * r^2
                V_sphere = 4/3 * pi *r^3;% Volume of that sphere = 4/3 * pi * r^3
                obj.v = obj.V/V_sphere;%%%  reduced volume
                obj.H = H;
                %obj.KG = KG;
                obj.Ro = sqrt(obj.A/4/pi);                     % radius of sphere of area A
                obj.m  = sum(H(:).*obj.basis.w.*SSn)/obj.Ro;   % reduced total mean curvature
        end
        function obj = update(obj)
            if obj.needs_updating
                obj.x = zeros(length(obj.basis.u(:)), 1);
                obj.y = zeros(length(obj.basis.u(:)), 1);
                obj.z = zeros(length(obj.basis.u(:)), 1);
                
                xu = zeros(length(obj.basis.u(:)), 1);
                yu = zeros(length(obj.basis.u(:)), 1);
                zu = zeros(length(obj.basis.u(:)), 1);
                
                xv = zeros(length(obj.basis.u(:)), 1);
                yv = zeros(length(obj.basis.u(:)), 1);
                zv = zeros(length(obj.basis.u(:)), 1);
                
                xuu = zeros(length(obj.basis.u(:)), 1);
                yuu = zeros(length(obj.basis.u(:)), 1);
                zuu = zeros(length(obj.basis.u(:)), 1);
                
                xvv = zeros(length(obj.basis.u(:)), 1);
                yvv = zeros(length(obj.basis.u(:)), 1);
                zvv = zeros(length(obj.basis.u(:)), 1);
                
                xuv = zeros(length(obj.basis.u(:)), 1);
                yuv = zeros(length(obj.basis.u(:)), 1);
                zuv = zeros(length(obj.basis.u(:)), 1);
                
                for mix = 1:obj.M,
                    for nix = 1:obj.N,
                        %%% calculate the coordinates
                        obj.x = obj.x + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        obj.y = obj.y + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        obj.z = obj.z + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the first derivatives w.r.t. u
                        xu = xu + (obj.ax(nix)*obj.basis.ducosnu(:,nix) + obj.bx(nix)*obj.basis.dusinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        yu = yu + (obj.ay(nix)*obj.basis.ducosnu(:,nix) + obj.by(nix)*obj.basis.dusinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        zu = zu + (obj.az(nix)*obj.basis.ducosnu(:,nix) + obj.bz(nix)*obj.basis.dusinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the first derivatives w.r.t. v
                        xv = xv + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvsinmv(:,mix));
                        yv = yv + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvsinmv(:,mix));
                        zv = zv + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvsinmv(:,mix));
                        
                        %%% calculate the second derivatives w.r.t. u
                        xuu = xuu + (obj.ax(nix)*obj.basis.duducosnu(:,nix) + obj.bx(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cx(mix)*obj.basis.cosmv(:,mix) + obj.dx(mix)*obj.basis.sinmv(:,mix));
                        yuu = yuu + (obj.ay(nix)*obj.basis.duducosnu(:,nix) + obj.by(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cy(mix)*obj.basis.cosmv(:,mix) + obj.dy(mix)*obj.basis.sinmv(:,mix));
                        zuu = zuu + (obj.az(nix)*obj.basis.duducosnu(:,nix) + obj.bz(nix)*obj.basis.dudusinnu(:,nix)).*(obj.cz(mix)*obj.basis.cosmv(:,mix) + obj.dz(mix)*obj.basis.sinmv(:,mix));
                        
                        %%% calculate the second derivatives w.r.t. v
                        xvv = xvv + (obj.ax(nix)*obj.basis.cosnu(:,nix) + obj.bx(nix)*obj.basis.sinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvdvsinmv(:,mix));
                        yvv = yvv + (obj.ay(nix)*obj.basis.cosnu(:,nix) + obj.by(nix)*obj.basis.sinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvdvsinmv(:,mix));
                        zvv = zvv + (obj.az(nix)*obj.basis.cosnu(:,nix) + obj.bz(nix)*obj.basis.sinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvdvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvdvsinmv(:,mix));
                        
                        %%% calculate duv
                        xuv = xuv + (obj.ax(nix)*obj.basis.ducosnu(:,nix) + obj.bx(nix)*obj.basis.dusinnu(:,nix)).*(obj.cx(mix)*obj.basis.dvcosmv(:,mix) + obj.dx(mix)*obj.basis.dvsinmv(:,mix));
                        yuv = yuv + (obj.ay(nix)*obj.basis.ducosnu(:,nix) + obj.by(nix)*obj.basis.dusinnu(:,nix)).*(obj.cy(mix)*obj.basis.dvcosmv(:,mix) + obj.dy(mix)*obj.basis.dvsinmv(:,mix));
                        zuv = zuv + (obj.az(nix)*obj.basis.ducosnu(:,nix) + obj.bz(nix)*obj.basis.dusinnu(:,nix)).*(obj.cz(mix)*obj.basis.dvcosmv(:,mix) + obj.dz(mix)*obj.basis.dvsinmv(:,mix));
                        
                    end
                end
                obj.x = reshape(obj.x,size(obj.basis.u));
                obj.y = reshape(obj.y,size(obj.basis.u));
                obj.z = reshape(obj.z,size(obj.basis.u));
                
                xu = reshape(xu,size(obj.basis.u));
                yu = reshape(yu,size(obj.basis.u));
                zu = reshape(zu,size(obj.basis.u));
                
                xv = reshape(xv,size(obj.basis.u));
                yv = reshape(yv,size(obj.basis.u));
                zv = reshape(zv,size(obj.basis.u));
                
                xuu = reshape(xuu,size(obj.basis.u));
                yuu = reshape(yuu,size(obj.basis.u));
                zuu = reshape(zuu,size(obj.basis.u));
                
                xvv = reshape(xvv,size(obj.basis.u));
                yvv = reshape(yvv,size(obj.basis.u));
                zvv = reshape(zvv,size(obj.basis.u));
                
                xuv = reshape(xuv,size(obj.basis.u));
                yuv = reshape(yuv,size(obj.basis.u));
                zuv = reshape(zuv,size(obj.basis.u));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% Calculate first and second fundamental forms%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                X  =[obj.x(:) obj.y(:) obj.z(:)];
                Xt = [xu(:) yu(:) zu(:)];
                Xp = [xv(:) yv(:) zv(:)];
                Xpp = [xvv(:) yvv(:) zvv(:)];
                Xtp = [xuv(:) yuv(:) zuv(:)];
                Xtt = [xuu(:) yuu(:) zuu(:)];
                
                E = dot(Xt,Xt,2);
                F = dot(Xt,Xp,2);
                G = dot(Xp,Xp,2);
                SS = (cross(Xt,Xp,2));SSn = sqrt(E.*G-F.*F);
                n = SS./SSn(:,ones(1,3));
                L = dot(Xtt,n,2);
                M = dot(Xtp,n,2);
                N = dot(Xpp,n,2);
                %%%%%%%%%%%%%%%%%%%%%%%% update geometrical properties
                obj.V = abs(1/3.*sum(obj.basis.w.*(dot(X,n,2)).*SSn));
                obj.A = sum(obj.basis.w.*SSn);
                H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
                KG = (L.*N - M.*M)./(E.*G-F.*F);
                %% quantities of interest
                obj.h = 1./obj.A.*sum(H(:).*obj.basis.w.*SSn);                	% normalized total mean curvature
                obj.T = sum(sum(KG(:).*obj.basis.w.*SSn))/2/pi;            		% total Gaussian curvature (constant for topology) should be zero for torus
                %obj.S = (2.*obj.H.^2-obj.KG).^(0.5);                     		% curvedness -- uncomment as needed
                obj.Eb = 1/2*sum((2.*H(:)).^2.*obj.basis.w.*SSn)./4/pi/pi;      % bending energy relative to that of the Clifford torus
                
                r = sqrt(obj.A/4/pi);% Area of a sphere = 4 * pi * r^2
                V_sphere = 4/3 * pi *r^3;% Volume of that sphere = 4/3 * pi * r^3
                obj.v = obj.V/V_sphere;%%%  reduced volume
                obj.H = H;
                obj.KG = KG;
                obj.Ro = sqrt(obj.A/4/pi);                     % radius of sphere of area A
                obj.m  = sum(H(:).*obj.basis.w.*SSn)/obj.Ro;   % reduced total mean curvature
            end
        end
        function plot(obj)
            if obj.needs_updating, obj = update(obj);end
            %% plot the shape
            surf(obj.x,obj.y,obj.z, 'EdgeColor', 'none');
            axis equal;
            view(3);
            lighting gouraud;
            camlight;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
            axis off;
            
        end
        
        function obj = set_X_o(obj, X)
            obj.X_o = X;
            N = obj.N;
            M = obj.M;
            nc = length(X)/3;
            xvec = X(1:nc); yvec = X(nc+1:2*nc);zvec = X(2*nc+1:end);
            obj.ax = xvec(1:N);obj.bx = xvec(N+1:2*N);obj.cx = xvec(2*N+1:(2*N)+M);obj.dx = xvec(((2*N)+ M+1):end);
            obj.ay = yvec(1:N);obj.by = yvec(N+1:2*N);obj.cy = yvec(2*N+1:(2*N)+M);obj.dy = yvec(((2*N)+ M+1):end);
            obj.az = zvec(1:N);obj.bz = zvec(N+1:2*N);obj.cz = zvec(2*N+1:(2*N)+M);obj.dz = zvec(((2*N)+ M+1):end);
            obj.needs_updating = 1;
        end
        function X_o = get_X_o(obj)
            X_o = [obj.ax(:);obj.bx(:);obj.cx(:);obj.dx(:);...
                   obj.ay(:);obj.by(:);obj.cy(:);obj.dy(:);...
                   obj.az(:);obj.bz(:);obj.cz(:);obj.dz(:)];
        end
        %         function obj = set_center(obj, pos)
        %             %%% position is given as absolute xyz
        %             obj.xc(1)   = pos(1);
        %             obj.yc(1)   = pos(2);
        %             obj.zc(1)   = pos(3);
        %             obj.X_o     = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
        %         end
        %         function obj = scale(obj, s)
        %             if length(s) == 1, s = [s s s];end
        %             obj.xc(2:end) = obj.xc(2:end)*s(1);
        %             obj.yc(2:end) = obj.yc(2:end)*s(2);
        %             obj.zc(2:end) = obj.zc(2:end)*s(3);
        %             obj.X_o     = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
        %         end
        %         function obj = rotate_shps_around_self(obj,ang)
        %             %%% rotate the shapes described by the spherical harmonic coefficiencts S
        %             %%% by the Euler angles a b g (radians) around each shape's center
        %             %%% rotation conventions are y-z-y
        %             %c = c/sqrt(1/4/pi);
        %             %% generate rotation matrix
        %             a = ang(1);b = ang(2);g = ang(3);
        %             Rg = rot_mx(g,2);
        %             Rb = rot_mx(b,1);
        %             Ra = rot_mx(a,2);
        %             R = Ra*Rb*Rg;
        %             %% loop over the shapes and rotate them around the center
        %             c = [obj.xc(1) obj.yc(1) obj.zc(1)];
        %             obj.xc(1) = 0;obj.yc(1) = 0;obj.zc(1) = 0;
        %             C = [obj.xc(:)';obj.yc(:)';obj.zc(:)'];
        %             cr = R*C;
        %             cr = cr';
        %             tx = [cr(:,1)];tx(1) = tx(1) + c(1);
        %             ty = [cr(:,2)];ty(1) = ty(1) + c(2);
        %             tz = [cr(:,3)];tz(1) = tz(1) + c(3);
        %             obj.X_o = [tx; ty; tz ];
        %         end
        %         function obj = mesh2shp(obj, m, L_max)
        %             if nargin<3,L_max = obj.L_max;disp(['Using L_max: ' num2str(L_max)]);end
        %             m = m.map2sphere;
        %             obj = obj.shp_analysis(m.X, m.t, m.p, L_max);
        %         end
        %         function obj = shp_analysis(obj, X, t, p, L_max)
        %             %%% The expansion of the three functions X = [x(t,p), y(t,p) z(t,p)] on the sphere.
        %             %%% [xyz] and the vectors are defined on the sphere at the spherical
        %             %%% coordinates t and p.
        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             [L, K ] = shp_surface.indices_gen(1:(L_max + 1)^2);
        %             M = length(L); %% number of functions in expansion
        %             N = length(X(:,1)); %% number of data points
        %             A  = zeros(N, M, 'double');
        %             for S = 1:length(L),
        %                 A(:,S) = obj.basis.ylk_bosh(L(S),K(S),p(:)',t(:)')';
        %             end;
        %             [U, S, V] = svd(A, 'econ');invS = 1./(S);invS(invS==inf) = 0;
        %             [clks] = (V*invS) * (U'*X);
        %             obj.X_o = double(reshape(clks,1,M*3));
        %         end
    end
    %     methods
    %         function obj = set.A(obj,A)
    %             obj.A = A;
    %         end
    %         function obj = set.V(obj,V)
    %             obj.V = V;
    %         end
    %         function obj = set.H(obj,H)
    %             obj.H = H;
    %         end
    %         function obj = set.KG(obj, KG)
    %             obj.KG = KG;
    %         end
    %         function obj = set.T(obj, T)
    %             obj.T = T;
    %         end
    %         function obj = set.h(obj, h)
    %             obj.h = h;
    %         end
    %         function obj = set.S(obj, S)
    %             obj.S = S;
    %         end
    %         function obj = set.Eb(obj, Eb)
    %             obj.Eb = Eb;
    %         end
    %     end
end















