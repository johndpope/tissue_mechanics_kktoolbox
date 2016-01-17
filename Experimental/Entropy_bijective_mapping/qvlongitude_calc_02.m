function [p, A, b, dl] = qlongitude_calc_02(X,t, A, F, L, ixN, ixS, fm, Ld, Li)
% Calculate longitude from diffusion as in Brechbuehler 1995
% Example:
%       p = longitude_calc(x, y, z, t, A, L, ixN, ixS)
% INPUT:
%       t: vector of latitude temperatures from lattitude_calc.m
%       A: the matrix A as obtained from lattitude_calc.m 
%       L: Cell array of link arrays
%       ixN and ixS: Identify the northpole and the south pole vertices
% OUTPUT:
%       p: phi value associated with each vertex
%       A and b: matrix A and b as calculated in the reference above
% Author:
%       Khaled Khairy  September 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the matrix A as in the pseudo-code given on page 158
% For the north pole
links = Ld{ixN};     % obtain the array of indices that the North pole is liked to
for ix = 1:length(links),       % loop over links
    A(links(ix),links(ix)) = A(links(ix),links(ix))-1;      % cut the link with the North pole
end
% For the south pole
links = Ld{ixS};     % obtain the array of indices that the south pole is liked to
for ix = 1:length(links),       % loop over links
    A(links(ix),links(ix)) = A(links(ix),links(ix))-1;      % cut the link with the south pole
end
A(1,1) = A(1,1) + 2;    % additional condition can be added to any row
b = sparse(length(Ld),1);% length L is the number of vertices
previous = ixN;
nbrs = Ld{ixN};
here = nbrs(1);       % any neighbor of the north pole
counter = 0;
maximum = 0;
dl = [ixN here];
%figure;u = [];
dt = [];
while(here~=ixS)
    %disp(counter);
    counter = counter + 1;if counter>length(L),error('could not determine p');end
    dt(counter) = here;
    nbrs = Ld{here}; % get the direct neighbors of here (array of neighbor indices)
    for ix = 1:length(nbrs),    % loop over neighbors of here
        if t(nbrs(ix)) > maximum,
            maximum = t(nbrs(ix));
            nextpos = ix;
        end
        if nbrs(ix) == previous
            prevpos = ix;
        end
    end
    % at this point both nextpos and prevpos are correct (I checked)
    % the problem is that we don't know who of the neighbors 
    % is "west" AND directly connected along
    % the path formed by previous->here->nbrs(nextpos),
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% debugging: plot the faces involved in the dateline
% %     u = unique([fm{here};fm{previous}]);  
% %     dfig(1);cla;patch(  'vertices',X,'faces',F(u,:),'FaceColor','interp',...
% %         'FaceVertexCData',b,'CDataMapping','scaled', 'EdgeColor','k');axis equal;
% %     pixs = F(u,:);pixs = unique(pixs);
% %     d = 0.1;
% %     for ix = 1:numel(pixs)
% %         str = sprintf(' %d',pixs(ix));
% %         text(X(pixs(ix),1)+d, X(pixs(ix),2)+d, X(pixs(ix),3)+d,str);
% %     end
% %     hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if previous==ixN,
        ihnbrs = Li{here};   % these are the indirect neighbors of "here"
        % there are four of them for a quadrangulation
        dpnbrs = Ld{previous};   % these are the direct neighbors of previous
        % there are four of them for a quadrangulation
        cnd = intersect(ihnbrs,dpnbrs);% determine the two vertices connected to "previous" directly and "here"
        % indirectly
        cc = cnd(1);
    end
    [W, cc]= get_cyclic_indices(cc,here,previous, nbrs(nextpos), fm, L, Ld, Li, F);
    disp([here W(:)']);
    b(W) = b(W) + 2 * pi; % 
    b(here) = b(here)- numel(W) * 2*pi;

    previous = here;
    here = nbrs(nextpos);
    dl = [dl here];
    % % % %     pl = kk_plot3([X(previous,:);X(here,:)]); set(pl,'LineWidth',5);hold on;
end

% %  solve the linear system
% [L1,U1] = luinc(A,1e-3);  %#ok<*REMFF1>
s.droptol = 1e-3;s.type = 'ilutp';[L1,U1] = ilu(A,s);
[p, flag1, relres1, iter1, resvec1] = bicg(A,b, 1e-6, 1000, L1, U1);

%% 
function [W,cc] = get_cyclic_indices(cc,here,previous, next, fm, L, Ld, Li, F)
    % W holds the indices as we go cyclically from previous to next.
    % the direction is determined based on cc.
    % Note: cc is connected directly to previous always and therefore sets 
    % a direction.
    W = [];
    p = L{here};    % these are the indices of all points around here
                    % in principle all we need is to sort those cyclically
                    % from previous to next and only consider direct
                    % connections
    % knowing cc, we know the shared face and can find the direct link to
    % here
        ccface = intersect(fm{here}, fm{cc});
        ix = F(ccface,:);
        ix(ix==here) = [];
        ix(ix==previous)=[];
        ix(ix==cc) = [];
        if ix==next,ccreturn = previous;end
    while(ix~=next)
        W = [W ix];
        %%% determine the new ccface
        newccface = intersect(fm{here}, fm{ix});    % gives two candidate faces (the ones that share the edge between here and ix
        newccface(newccface==ccface) = [];      % only leaves the new face
        ccface = newccface;
        %%% determine the new cc
        ccreturn = ix;
        cc = F(ccface,:);
        cc(cc==here) = [];
        cc(cc==ix)=[];
        cc(cc==intersect(Ld{here},cc)) = [];
        %%% determine index of direct connection
        ixold = ix;
        ix = F(ccface,:);
        ix(ix==here) = [];
        ix(ix==cc) = [];
        ix(ix==ixold) = [];
    end

    cc = ccreturn;




















