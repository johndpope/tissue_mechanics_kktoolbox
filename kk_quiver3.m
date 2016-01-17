function kk_quiver3(X1, X2, X3,D1,D2,D3, step)
if nargin == 6,step = 1;end
X1 = X1(1:step:length(X1));
X2 = X2(1:step:length(X2));
X3 = X3(1:step:length(X3));
D1 = D1(1:step:length(D1));
D2 = D2(1:step:length(D2));
D3 = D3(1:step:length(D3));

% h = quiver3(X1, X2, X3,D1,D2,D3,'b');
% adjust_quiver_arrowhead_size(h, 1.5);   % Makes all arrowheads 50% bigger.

c = mat2gray([D1/max(D1),D2/max(D2),D3/max(D3)]);
c = mat2gray(1./c);
for ix = 1:numel(X1)
% h = quiver3(X1(ix), X2(ix), X3(ix),D1(ix),D2(ix),D3(ix), 5,'Color',c(ix,:));
% h = quiver3(X1(ix), X2(ix), X3(ix),D1(ix),D2(ix),D3(ix), 9,'b');
% adjust_quiver_arrowhead_size(h, 12);   % scales arrowheads.
h = quiver3(X1(ix), X2(ix), X3(ix),D1(ix),D2(ix),D3(ix), 100 * c(ix),'w');
adjust_quiver_arrowhead_size(h, 10);   % scales arrowheads.
set(h,'linewidth',2);
end