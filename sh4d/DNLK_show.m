function DNLK_show(DNLK)
%% shows (and analyses) the coefficients categorized according to their component
%% amplitudes and orders in frequency
figure;
D = DNLK.D;
N = DNLK.N;
L = DNLK.L;
K = DNLK.K;

dim = length(D);
N = N(1:dim);L = L(1:dim);K = K(1:dim);


%%%%% show the individual Coefficient values in the order they come with
[ds, IX] = sort(D.^2);
subplot(4,1,1);plot(D.^2);

nshow = 1;

for ix = 0:nshow-1,
    str = sprintf('N%d\tL%d\tK%d\tD = %.5f',N(IX(end-ix)), L(IX(end-ix)), K(IX(end-ix)), D(IX(end-ix)));
    
    disp(str);
%     xpos = IX(end-ix);ypos = D(IX(end-ix)).^2; text(xpos,ypos,str);
end

%% generate boxplot for N
icounter = 0;
jcounter = 0;
N_vec = zeros(length(min(N):max(N)), length(N), 'double');
lab_cell = {};
for ix = min(N):max(N), % loop over possible N values
    icounter = icounter + 1;
    for jx = 1:length(N),
        if N(jx) == ix,
            jcounter = jcounter + 1;
            N_vec(icounter, jcounter) = D(jx);
        end
    end
    lab_cell{icounter} = num2str(ix);
end

subplot(4,1,2);boxplot(N_vec', 'labels', lab_cell);drawnow

%% generate boxplot for L
icounter = 0;
jcounter = 0;
L_vec = zeros(length(min(L):max(L)), length(L), 'double');
lab_cell = {};
for ix = min(L):max(L), % loop over possible L values
    icounter = icounter + 1;
    for jx = 1:length(L),
        if L(jx) == ix,
            jcounter = jcounter + 1;
            L_vec(icounter, jcounter) = D(jx);
        end
    end
    lab_cell{icounter} = num2str(ix);
end
subplot(4,1,3);boxplot(L_vec', 'notch', 'on', 'labels', lab_cell);drawnow


%% generate boxplot for K
icounter = 0;
jcounter = 0;
K_vec = zeros(length(min(K):max(K)), length(K), 'double');
lab_cell = {};
for ix = min(K):max(K), % loop over possible L values
    icounter = icounter + 1;
    for jx = 1:length(K),
        if L(jx) == ix,
            jcounter = jcounter + 1;
            K_vec(icounter, jcounter) = D(jx);
        end
    end
    lab_cell{icounter} = num2str(ix);
end
subplot(4,1,4);boxplot(K_vec', 'labels', lab_cell);drawnow






