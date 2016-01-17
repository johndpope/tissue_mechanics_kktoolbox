function [t, A, b] = qlatitude_calc(L, ixN, ixS)
% Calculate latitude from diffusion as in Brechbuehler 1995
% Example:
%        [t, A, b] = latitude_calc(L, ixN, ixS)
% INPUT:
%       L: Cell array of link arrays
%       ixN and ixS: Identify the northpole and the south pole vertices
% OUTPUT:
%       t: theta value associated with each vertex
%       A and b: the matrix A and vector b (as in Brechbuehler 1995)
% Author:
%       Khaled Khairy  September 2004
%       This function has been modified to work for triangulations also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up matrix A (as in the reference above on p.157 and the appendix)
A = sparse(length(L));     % creates a square zero matrix of size nv
b = sparse(length(L),1);
for iv = 1:length(L)      % loop over the cells inLF
    links = L{iv};  % the array that contains the indices of the vertices linked to iv
    A(iv,iv) = length(links);
    for K = 1:length(links)         % loop over the indices that are linked
        A(iv,links(K)) = -1;               % set all linked elements to -1
    end
end

A(ixN,:) = 0;A(ixN,ixN) = 1;A(ixS,:) = 0;A(ixS,ixS) = 1;
b(ixS) = pi;% Setup the vector b
%[L1,U1] = luinc(A,1e-6);
s.droptol = 1e-6;s.type = 'ilutp';[L1,U1] = ilu(A,s);
[t] = bicg(A,b, 1e-6, 100, L1, U1);
