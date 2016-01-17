function plot_time_sequence(X_o_res_vec,s_master,phys,pos, Y_LK_vm, az, el)

if nargin <= 5, az = -40;el = 0;end
zbnd = phys.zbnd;
close all;dfig;  
view(az,el);axis equal;axis on;
for ix = 1:size(X_o_res_vec,1)-1
    disp(pos(ix,:));
      xlim([-4 4]);ylim([-4 4]);zlim([-1.2 3]);
      X_slave = X_o_res_vec(ix,:)';
    [xclks yclks zclks] = get_xyz_clks(X_slave);
    X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];%
    
    cla;view(az,el)
    patch('vertices',X,'faces',phys.C_slave,'FaceColor','red', 'FaceAlpha', 0.8);
    hold on;
    
    
    
    s_master = position(s_master,pos(ix,:));
    phys.X_o_vm = s_master.X_o;
    [xc yc zc] = get_xyz_clks(phys.X_o_vm);
    phys.X_vm = Y_LK_vm(:,1:length(xc))* [xc(:) yc(:) zc(:)];% get vertices coordinates for vitteline membrane
    
    patch('vertices',phys.X_vm,'faces',phys.C_vm,'FaceColor','green', 'FaceAlpha', 1, 'EdgeColor','none');
    hold off
    lighting phong;camlight;
    drawnow
    
    [res_self] = tri_tri_self_intersect(X,phys.C_slave, phys.TP_self);
    [res_vm] = tri_tri_intersect(X,phys.C_slave, phys.X_vm, phys.C_vm, phys.TP_vm);
    if any(X(:,3)<zbnd), zviolation = 1e0;else zviolation = 0;end
    disp([res_self res_vm zviolation]);
%           pause
end
rotate3d on;