function ix = get_index(X,v)
dist = sqrt(sum([(X(:,1)-v(1)).^2 + (X(:,2)-v(2)).^2 + (X(:,3)-v(3)).^2],2));
[B,IX] = sort(dist);
ix = IX(1);
% % X1 = X;
% % if nargin<3, d = 0.04;end
% % ix = find(X(:,1)<val(1)+d & X(:,1)>val(1)-d);
% % if numel(ix)==1, return;
% % else
% %     X1 = X1(ix,:);
% %     ix = find(X(:,2)<val(2)+d & X(:,2)>val(2)-d);
% % end

