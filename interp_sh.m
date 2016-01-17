function [CLKS] = interp_sh(sclks,eclks,n)
%%% interpolates linearly between shape coefficients
%%% this generates a morph effect.
%%%     CLKS = interp_sh(sclks,eclks,n)
%%% and interpolating the coefficient values
%%% Returns CLKS, an n x length(sclks) matrix

sclks = sclks(:)';eclks = eclks(:)';
%% make sure both are same length as sclks
nc = round(length(sclks)/3);sxclks = sclks(1:nc);syclks = sclks(nc+1:2*nc); szclks = sclks(2*nc+1:3*nc);
L_max = sqrt(length(sclks)/3)-1;
nc = round(length(eclks)/3);exclks = eclks(1:nc);eyclks = eclks(nc+1:2*nc); ezclks = eclks(2*nc+1:3*nc);
nL_max = L_max;nc = (nL_max +1)^2;
nxclks = zeros(1,nc);nyclks = nxclks ; nzclks = nxclks;
if length(nxclks)>length(exclks), 
    nxclks(1:length(exclks)) = exclks;nyclks(1:length(exclks)) = eyclks;nzclks(1:length(exclks)) = ezclks;
else, 
    nxclks(:) = exclks(1:length(nxclks));nyclks(:) = eyclks(1:length(nyclks)); nzclks(:) = ezclks(1:length(nzclks));
end
eclks = [nxclks(:);nyclks(:);nzclks(:)]';
if length(sclks)~=length(eclks),error('starting and end shape must have same Lmax');end;
%%% interpolate between sclks and eclks to generate the required sequence
CLKS = [];
for ix = 1:length(sclks),  % loop over the coefficients
    clkinter = sclks(ix):sign(eclks(ix)-sclks(ix))*abs(eclks(ix)-sclks(ix))/n:eclks(ix);
    if sclks(ix)==eclks(ix),clkinter = ones(1,n+1)*sclks(ix);end;
    CLKS(:, ix) = clkinter(:);
end
CLKS = CLKS(2:end,:);