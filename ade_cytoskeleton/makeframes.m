function [v, da_o, w] = makeframes(L_max, gdim, p, t, L, K, phasename, pflag)
str = sprintf('load %s;',phasename);eval(str);
% if  sum(v_o_rec-v_o_rec(1))==0,
%     %% then we have a constant v_o rec
%     subplot(2,1,1);plot(da_o_rec, wb_rec,'b-');
%     xlabel('da_o'); ylabel('bending energy');
%     subplot(2,1,2);
% end
% if sum(da_o_rec-da_o_rec(1))==0,
%     %% then we have a constant da_o rec
%     subplot(2,1,1);plot(v_o_rec,wb_rec,'b-');
%     xlabel('v_o'); ylabel('bending energy');
%     subplot(2,1,2);
% end
da_o = da_o_rec;
w = wb_rec;
v = v_o_rec;
for ix = 1:length(wb_rec),
    %% check if the we have a constant v_o or a constant da_o
    
    xp  = x_vec(ix,:);
    str = sprintf('Bending energy = %.6f   v_o = %.6f  da_o = %.6f',...
        wb_rec(ix), v_o_rec(ix), da_o_rec(ix));
    h = length(xp);xp = xp(ones(gdim^2,1),:);xp = reshape(xp,gdim, gdim, h); % expensive
    plot_spharm(gdim, p, t, L, K, xp, str);
    view(3);
    xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
    light;lighting phong; axis('square');
    camlight headlight;
    grid off
    drawnow
    if pflag
        step = ix;
        if step<10,
            str = sprintf('L%d_gdim%d_000%d',L_max, gdim, step);
        elseif step<100
            str = sprintf('L%d_gdim%d_00%d',L_max, gdim, step);
        elseif step<1000
            str = sprintf('L%d_gdim%d_0%d',L_max, gdim, step);
        else
            str = sprintf('L%d_gdim%d_%d',L_max, gdim, step);
        end
        str = sprintf('print -dtiff -r600 %s%s.tif;',phasename,str);eval(str);
    else
%         pause(.5);
    end
end