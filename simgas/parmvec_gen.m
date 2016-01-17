function [parmvec phys] = parmvec_gen(X_o,phys)
% This function is under construction!! it should be a complete system to
% precisely define the fitting parameters and which parameters remain fixed


L_max = shp_surface.get_L_max(X_o);


%%% define the indices to be fixed. Use [] to fit all parameters
if isempty(phys.del_indx)
    phys.del_indx = [1 (L_max + 1)^2+1 2*(L_max + 1)^2 + 1];    % fix x y z
%phys.del_indx = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys.bin_parms = ones(size(X_o));
phys.bin_parms(phys.del_indx) = 0;
phys.del_parms = X_o(phys.del_indx);
%specify parmvec
parmvec = X_o;
parmvec(phys.del_indx) = [];
% end parmvec specification