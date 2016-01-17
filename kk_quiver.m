function kk_quiver(X2, X3,D2,D3, step)
if nargin == 4,step = 1;end

X2 = X2(1:step:length(X2));
X3 = X3(1:step:length(X3));
D2 = D2(1:step:length(D2));
D3 = D3(1:step:length(D3));

c = mat2gray([D2/max(D2),D3/max(D3)]);

for ix = 1:numel(X2)
 h = quiver(X2(ix), X3(ix),D2(ix),D3(ix), 5,'k');
adjust_quiver_arrowhead_size(h, 6);   % scales arrowheads.
set(h,'linewidth',1.5);
end