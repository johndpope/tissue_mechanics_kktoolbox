function [T] = clean_T(T, A_vec, V_vec, wb_vec)
% discard excessively large or strange structures by only keeping those
% that "make sense"


L_max = get_L_max(T(1,:));
wb_thresh = 8.0;
V_max = inf;
V_min = 0;
indx = find(wb_vec>wb_thresh);
for ix = 1:length(indx)
    X_o = T(indx(ix),:);
    X_o = tr(X_o,1);
    T(indx(ix),:) = tr(X_o,L_max);
end


%% % delete objects with the wrong volume
A_max = median(A_vec) * 3 ;  % upper limit on area
A_min = 0;
indx = [];
for ix = 1:length(A_vec)
    if ~(V_vec(ix)<A_max && A_vec(ix)>A_min),
        indx = [indx ix];
    end
end
T(indx,:) = [];
wb_vec(indx) = [];
V_vec(indx) = [];
A_vec(indx) = [];

disp(['Discarded ' num2str(length(indx)) ' objects based on area']);


%%

% % fac = 6;
% %
% % % med = median(A_vec);
% % % delix = find(A_vec<fac*med);
% %
% % delix = wb_vec<fac;
% %
% % V_vec_clean = V_vec(delix);
% % A_vec_clean = A_vec(delix);
% % wb_vec_clean = wb_vec(delix);
% % T_clean = T(delix, :);
