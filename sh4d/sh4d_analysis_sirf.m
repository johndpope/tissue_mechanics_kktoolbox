function [DNLK f] = sh4d_analysis_sirf(I,N_max,L_max)
% SIRF: Simple Iterative Residual Fitting
% Series expansion of the 4-D image (3D volume + intensity) in 4D spherical
% harmonic functions. Basis functions are generated by sh4d_basis_gen.
% Returns the coefficients from which the complete image may be
% reconstructed, in a low-pass filtered sense.

verbose = 0;
gran = 5;
dim = size(I,1);
tau = 100;%size(I, 1);
pad_vec = [0 0 0];%([dim dim dim]);
I = padarray(I,pad_vec);

[av_i, av_j, av_k] = center_of_mass(I);
av_i = round(av_i);av_j = round(av_j);av_k = round(av_k);
%% generate theta, phi and r from the image
[ix iy iz] = ind2sub(size(I),1:length(I(:)));
x = ix -av_i;y = iy-av_j;z = iz-av_k;  % move coordinate system to center of mass
[t p r] = kk_cart2sph(x,y,z);
indx = sub2ind(size(I),ix, iy, iz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the first set of required basis vectors

[N, L, K ] = sh4d_indices_gen(N_max,L_max);

% % % % %%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%
% % % % N_ex = [];  % N values to be excluded
% % % % for ix = 1:length(N_ex),
% % % %     L(N==N_ex(ix)) = [];K(N==N_ex(ix)) = [];N(N==N_ex(ix))=[];
% % % % end
% % % % %%%%%%% Exclude basis vectors for L here
% % % % L_ex = [1:6];  % L values to be excluded
% % % % for ix = 1:length(L_ex),
% % % %     N(L==L_ex(ix)) = [];K(L==L_ex(ix)) = [];L(L==L_ex(ix))=[];
% % % % end
% % % % %%%%%%% Exclude basis vectors for K here
% % % % K_ex = [1:6];%[[-max(L):-1] [1:5]];  % K values to be excluded
% % % % for ix = 1:length(K_ex),
% % % %     N(K==K_ex(ix)) = [];L(K==K_ex(ix)) = [];K(K==K_ex(ix))=[];
% % % % end
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% % % % %%% sort the indices to process low frequencies first
%  N(K<0) = [];L(K<0) = [];K(K<0) = [];
bv = [N L K];
[B,IX] = sort(bv(:,3).^2); bv = (bv(IX,:));  % sort according to K
N = bv(:,1);L = bv(:,2);K = bv(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Total number of basis functions: ' num2str(length(N))]);
% s = length(N);
maxDim = gran;
if length(N)> maxDim, s = maxDim;
else
    s = length(N);
end

n = length(r(:)); %% number of data points
str = sprintf('SVD_gdim%d_N%d_L%d_s%d.mat',size(I,1),N_max, L_max, s);
if exist(str)==2,
    load(str);
    if verbose, disp('Loading precalculated basis matrix decomposition');end;
else
    warning off MATLAB:divideByZero;
    A  = zeros(n, s, 'single');
    h = waitbar(0,'Generating basis matrix ....');
    for S = 1:s,
        A(:, S) = sh4d_basis_gen(N(S),L(S),K(S),p(:)',t(:)',r(:)', tau)';
        waitbar(S/s);
    end
    close(h);
    drawnow;
    disp('Solving initial linear system...');
    [U, S, V] = svd(A, 0);
    invS = 1./(S);invS(invS==inf) = 0;
%     save(str,'U', 'V', 'invS', 'A');
    disp('Done !');
end

%%% use the svd
DNLK.D= (V*invS) * (U'*(I(:)));
disp('Initial LS calculation ..... done');

%% Calculate the residual
f = A*DNLK.D;
% kk_montage(mat2gray(reshape(f,size(I))));drawnow;
clear A U S V invS;
disp('Cleanup ..... done');

% % %%% just in case you want to look at it at this stage
% % im = I;
% % im(indx) = f;
% % %% crop the arrays
% % im = im(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
% % I = I(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
% % %% display
% % im = mat2gray(im);
% % kk_montage([mat2gray([I]) (im)]);
% % impixelinfo;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iteratively refine residual --- do sirf here
residual = I(:)-f;
h = waitbar(0,'SIRF iterations in progress ....');
counter = 0;wl = length(s+1:length(L));
g = gran;%(L_max+1)^2;        % granularity for sirf
A = zeros(n,g, 'single');
res_vec = [];
dfig;
for S = s+1:g:length(L),
    disp([N(S) L(S) K(S)]);
    acounter = 0;

    for iS = S:S+g-1
        if length(N)>=iS,
            acounter = acounter + 1;
            
            A(:,acounter) = sh4d_basis_gen(N(iS),L(iS),K(iS),p(:)',t(:)',r(:)', tau)'; % get basis (from database in future)
            %         A = sh4d_basis_gen(N(S),L(S),K(S),p(:)',t(:)',r(:)', tau)'; % getbasis (from database in future)
        end
    end
%   A = A(:,acounter);
    D = A\residual(:);
    %[U, SS, V] = svd(A,0); warning off; invS = 1./(SS);invS(invS==inf) = 0;warning on; D    = (V*invS) * (U'*(residual(:)));
    r_vec = A*D(:);
    r_vec = r_vec(:);
    f = f + r_vec;
    residual = residual-r_vec;
    res_vec = [res_vec (norm(I(:)-f(:)))];
    DNLK.D = [DNLK.D;D(:)];
    counter = counter+1;
    waitbar(counter/wl);
    temp = reshape(f,size(I));
    warning off;
    kk_montage([mat2gray(I(pad_vec(1)+1:end-pad_vec(1), pad_vec(2)+1:end-pad_vec(2), pad_vec(3)+1:end-pad_vec(3)))...
        mat2gray(temp(pad_vec(1)+1:end-pad_vec(1), pad_vec(2)+1:end-pad_vec(2), pad_vec(3)+1:end-pad_vec(3)))]);warning on;
    title('SIRF Analysis');
    
    drawnow;
end
close(h);


disp('...done');
f = reshape(f,size(I));
f = f(pad_vec(1)+1:end-pad_vec(1), pad_vec(2)+1:end-pad_vec(2), pad_vec(3)+1:end-pad_vec(3));
I = I(pad_vec(1)+1:end-pad_vec(1), pad_vec(2)+1:end-pad_vec(2), pad_vec(3)+1:end-pad_vec(3));


DNLK.D = DNLK.D(1:length(N));
DNLK.N = N;
DNLK.L = L;
DNLK.K = K;
DNLK.method = 'sirf';
DNLK.tau = tau;
DNLK.old = 0;
DNLK.im = mat2gray(f);

% kk_montage([mat2gray(I) mat2gray(f)]);
% impixelinfo;

disp(DNLK);DNLK_show(DNLK);

%% display

if counter, figure;plot(res_vec);end
%% Plot the volume based on the DNLKs only
% im = sh4d_plot(DNLK, N_max, L_max, 32);

% % im = im(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
% % I = I(pad_vec(1):end-pad_vec(1), pad_vec(2):end-pad_vec(2), pad_vec(3):end-pad_vec(3));
% % %% display
% % kk_montage(mat2gray([I im]));
% % impixelinfo;


