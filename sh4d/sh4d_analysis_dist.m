function [DNLK, im,A] = sh4d_analysis_dist(I,N_max,L_max)
% Series expansion of the 4-D image (3D volume + intensity) in 4D spherical
% harmonic functions. Basis functions are generated by sh4d_basis_gen.
% Returns the coefficients from which the complete image may be
% reconstructed, in a low-pass filtered sense.
% Note on memory requiremens: a 64^3 column requires 1MB memory,
% To calculate the memory requirements for the whole basis matrix do this:
% (xdim)x (ydim)x(zdim)x 4 x (no.of basis vetors)
% where number of basis vectors is :
%
% [N, L, K ] = sh4d_indices_gen(Nmax,Lmax); disp(length(L))
%
% for 16 GB total RAM a basis of 15 000 should be maximum (will try this).
% Nmax = 15, Lmax = 20 seems a good choice, but always keep s < 15000.

verbose = 1;
dim = round(size(I,1)/2);
pad_vec = [0 0 0];%([dim dim dim]);
I = padarray(I,pad_vec);

tau = size(I, 1);
% [av_i, av_j, av_k] = center_of_mass(I);av_i = round(av_i);av_j = round(av_j);av_k = round(av_k);
av_i = round(size(I,1)/2);av_j = round(size(I,2)/2);av_k = round(size(I,3)/2);

%% generate theta, phi and r from the image
[ix iy iz] = ind2sub(size(I),1:length(I(:)));
x = ix -av_i;y = iy-av_j;z = iz-av_k;  % move coordinate system to center of mass
[t p r] = kk_cart2sph(x,y,z);


%% generate the required basis vectors
%%% check whether we have the matrix A or svd components stored on disk
[N, L, K ] = sh4d_indices_gen(N_max,L_max);
s = length(L)-mod(length(L),4);s = round(s/4);
% if length(L)<s, s = length(L);end;

str = sprintf('SVD_gdim%d_N%d_L%d_%d.mat',size(I,1),N_max, L_max,s);
if exist(str)==2,
    load(str);
    if verbose, disp('Loading precalculated basis matrix decomposition');end;
else
    
    warning off MATLAB:divideByZero;
    %[N, L, K ] = sh4d_indices_gen(N_max,L_max);
    m = length(L); %% number of functions in expansion
    disp(m);
    n = length(r(:)); %% number of data points

    %%%% Calculate the basis matrix using nproc processors
        disp('Using 4 processors for building basis matrix');
        %%% generate input arguments for individual tasks
        
        vec = 1:s;    inp1{1} = N(vec);inp1{2} = L(vec);inp1{3} = K(vec);inp1{4} = p(:)';inp1{5} = t(:)';inp1{6} = r(:)';inp1{7} = tau;
        vec = s+1:2*s;inp2{1} = N(vec);inp2{2} = L(vec);inp2{3} = K(vec);inp2{4} = p(:)';inp2{5} = t(:)';inp2{6} = r(:)';inp2{7} = tau;
        vec = 2*s+1:3*s;inp3{1} = N(vec);inp3{2} = L(vec);inp3{3} = K(vec);inp3{4} = p(:)';inp3{5} = t(:)';inp3{6} = r(:)';inp3{7} = tau;
        vec = 3*s+1:4*s;inp4{1} = N(vec);inp4{2} = L(vec);inp4{3} = K(vec);inp4{4} = p(:)';inp4{5} = t(:)';inp4{6} = r(:)';inp4{7} = tau;

        sched = findResource('scheduler', 'type', 'local');
        j = createJob(sched);
                     
        warning off;
        createTask(j, @sh4d_basis_gen_vec, 1, {inp1});
        createTask(j, @sh4d_basis_gen_vec, 1, {inp2});
        createTask(j, @sh4d_basis_gen_vec, 1, {inp3});
        createTask(j, @sh4d_basis_gen_vec, 1, {inp4});
%         createTask(j, @sh4d_basis_gen_vec, 1, {inp1 inp2 inp3 inp4});
        warning on;
        
        %%% before submission test all inputs individually
% %         H = sh4d_basis_gen_vec(inp4{1}, inp4{2}, inp4{3}, inp4{4}, inp4{5}, inp4{6}, inp4{7});
% %         H = sh4d_basis_gen_vec(inp3{1}, inp3{2}, inp3{3}, inp3{4}, inp3{5}, inp3{6}, inp3{7});
% %         H = sh4d_basis_gen_vec(inp2{1}, inp2{2}, inp2{3}, inp2{4}, inp2{5}, inp2{6}, inp2{7});
% %         H = sh4d_basis_gen_vec(inp1{1}, inp1{2}, inp1{3}, inp1{4}, inp1{5}, inp1{6}, inp1{7});
  
        
        submit(j);
        waitForState(j, 'finished');
        
        %%%% display errors if any
        errmsgs = get(j.Tasks, {'ErrorMessage'});
        nonempty = ~cellfun(@isempty, errmsgs);
        celldisp(errmsgs(nonempty));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        results = getAllOutputArguments(j)';disp(results);
        if isempty(results{end}),
            error('Results dimensions inconsistency');
        end
        A = cell2mat(results);
        clear results;
        disp('basis matrix size');
        disp(size(A));
        if size(A)~=[length(p) 4*s],
            error('Basis matrix dimensions inconsistency');
        end
% %     else
% %     error('mod(s,nproc) is not equal to zero');
% %     end
    disp('Solving LS system....');
    warning off;
    [U, S, V] = svd(A'*A, 0);
    invS = 1./(S);invS(invS==inf) = 0;warning on;
    save(str,'U', 'V', 'invS', 'A');
end;

%%% use the svd
D= (V*invS) * (U'*(A'*I(:)));

disp('Done !');
%% Plot the volume based on the DNLKs only
f = A*D;
indx = sub2ind(size(I),ix, iy, iz);
im = I;
im(indx) = f;
%% crop the arrays
if dim,
    im = im(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
    I = I(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
    %% display
end
im = mat2gray(im);
kk_montage([mat2gray([I]) (im)]);
impixelinfo;
DNLK.D = D;
DNLK.N = N;
DNLK.L = L;
DNLK.K = K;
DNLK.tau = tau;



















