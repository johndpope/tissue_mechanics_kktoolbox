function [T, A_vec, V_vec, wb_vec]  =...
    final_clean_T(T, A_vec, V_vec, wb_vec, A_max, V_max, E_max)
%%% we need to get rid of the spurious shapes in T by replacing the  the
%%% spurious shapes with spheres of volume similar to neighbors
if mod(size(T,2), 3)>0, T = T(:,1:end-1);end
L_max = get_L_max(T(1,:));
medVol = median(V_vec);
R = (medVol * 3/4/pi)^(1/3);
A = 4 * pi * R^2;
X_o = 4*over_Basis(tr([0 0 0 1 0 1 0 0 0 0 1 0], L_max));
count = 0;
%% get rid of shapes that don't fall within the constraints
for ix  = 1: size(T,1),         %loop over shapes
    if A_vec(ix)>A_max || V_vec(ix)>V_max || wb_vec(ix)>E_max,
        [xc1 yc1 zc1] = get_xyz_clks(X_o);
        [xc2 yc2 zc2] = get_xyz_clks(T(ix,:));
        xc1(1) = xc2(1);yc1(1) = yc2(1);zc1(1) = zc2(1);
        T(ix,:) = [xc1(:)' yc1(:)' zc1(:)'];
        A_vec(ix) = A;
        V_vec(ix) = medVol;
        wb_vec(ix) = 1;
        count = count + 1;
    end
end
disp(['Replaced ' num2str(count) ' of ' num2str(size(T,1))]);
