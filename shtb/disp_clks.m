function disp_clks(clks, plot_flag)
%%% needs to be developed to distinguish different L orders

if nargin == 1, plot_flag = 0;end
[xc yc zc] = get_xyz_clks(clks);

if plot_flag
clf
subplot(3,1,1); plot(xc, '*');
subplot(3,1,2); plot(yc, '*');
subplot(3,1,3); plot(zc, '*');
end
L_max = get_L_max(clks);
last = 0;
for ix = 0:L_max
    disp(['L = ' num2str(ix)]);
    start = last+1; finish = (ix + 1).^2;
    vec = start:finish;
    last = finish;
    %vec
    disp([xc(vec) yc(vec) zc(vec)]);
end
