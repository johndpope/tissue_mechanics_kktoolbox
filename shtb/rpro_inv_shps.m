function [S] = rpro_inv_shps(S)
%%% SHP invariant description rotating parameter space and object space
if mod(size(S,2),3)>0, S = S(:,1:end-1);end
for ix = 1:size(S,1),   % loop over shapes
    disp(['Processing shape ' num2str(ix) ' of ' num2str(size(S,1))]);
    S(ix,:) = rpro_inv(S(ix,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
