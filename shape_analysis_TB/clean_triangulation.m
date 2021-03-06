function [X2, F2] = clean_triangulation(X, F, Fa)
%%% Cleans the triangulation from repetitive vertices
%%% and faces with zero area.
ind = 1:length(F);Fr = find(Fa==0);   % finds the faces that will be removed
ind = ind(find(Fa~=0));
% disp([find(Fa==0)]);        % these are the faces that will be removed
F2 = F(ind,:);

%%% now remove the duplicate vertices
nV = (length(F2)+4)/2;      %
%% Build a list of vertices to be removed


X2 = X;R = [];
for ix = 1:length(X2),
    xd = X2(ix,:);
    Xd = xd(ones(length(X2),1),:);
    aa  = X2==Xd;dup = sum(aa,2)==3;dup = find(dup);    % finds the rows that are identical
    if length(dup)>1,R = [R ix];disp([ix X(ix,1) X(ix,2) X(ix,3)]); end
end
%%



%% Manually determine the indices to keep and the ones to change
change = [15 13;54 13;1183 1058;3282 3227];
for ix = 1:size(change,1),
    F2(F2==change(ix,1)) = change(ix,2);
end

%% Now the indices that should be removed can be safely removed

[X2, F2]  = remove_unused_vertices(X2,F2,change(:,1));
% X2 = X;
% for ix = 1:length(X2)
%     xd = X2(ix,:);
%     Xd = xd(ones(length(X2),1),:);
%     aa  = X2==Xd;
%     dup = sum(aa,2)==3;    % finds the rows that are identical
%     dup = find(dup);
%     if isempty(dup), error('Vertex not found');end
%     if  length(dup)>1, 
%         [X2,F2] = remove_unused_vertices(X2,F2,dup(2:end));
%         disp(dup(2:end));
%     end
% end
