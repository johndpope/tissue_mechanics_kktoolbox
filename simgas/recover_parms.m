function [X_o] = recover_parms(parmvec,phys)
X_o = phys.bin_parms;
X_o(X_o==1) = parmvec;
X_o(phys.del_indx) = phys.del_parms;



% % phys.del_indx = [1 (L_max + 1)^2+1];
% % phys.bin_parms = ones(size(X_o));
% % phys.bin_parms(phys.del_indx) = 0;
% % phys.del_parms = X_o(phys.del_indx);
% % %specify parmvec
% % parmvec = X_o;
% % parmvec(phys.del_indx) = [];
% % % end parmvec specification