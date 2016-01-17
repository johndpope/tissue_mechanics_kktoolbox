function [Rx Ry Rz] = analyze_clks_difference(X_o1, X_o2, L_max, thresh)
% returns a filtered cell array. Each cell containts:
% label, L , K, val
plotflag = 0;
if numel(X_o1)~= numel(X_o2), error('vectors must be of same length');end

[xc1 yc1 zc1] = get_xyz_clks(X_o1);
[xc2 yc2 zc2] = get_xyz_clks(X_o2);

if nargin <3, L_max = get_L_max(X_o1);end
if nargin <4, thresh = 70;disp(['using threshold of ' num2str(thresh) ' percent']);end
if plotflag,
    dfig(1);clf;set(gcf,'Color','Black');
    dfig(2);clf;set(gcf,'Color','Black');
    dfig(3);clf;set(gcf,'Color','Black');
    plot_count = 1;
end
counter = 1;
Rx = {};
Ry = {};
Rz = {};
for l = 0:L_max
    lab = label_gen(l);
    lxc1 = [];
    lyc1 = [];
    lzc1 = [];
    lxc2 = [];
    lyc2 = [];
    lzc2 = [];
    for m = -l:l
        lxc1 = [lxc1 xc1(counter)];
        lyc1 = [lyc1 yc1(counter)];
        lzc1 = [lzc1 zc1(counter)];
        lxc2 = [lxc2 xc2(counter)];
        lyc2 = [lyc2 yc2(counter)];
        lzc2 = [lzc2 zc2(counter)];
        counter = counter + 1;
    end
    maxxc = max(abs([lxc1 lxc2]));
    maxyc = max(abs([lyc1 lyc2]));
    maxzc = max(abs([lzc1 lzc2]));
    
    maxc = max([maxxc maxyc maxzc]);
    disp(['L = ' num2str(l) ' maxcoeff = ' num2str(maxc)]);
    if(maxc>0)
        dlxc = abs(lxc1-lxc2)/maxc * 100;
        dlyc = abs(lyc1-lyc2)/maxc * 100;
        dlzc = abs(lzc1-lzc2)/maxc * 100;
    else
        dlxc = zeros(numel(lxc1));
        dlyc = zeros(numel(lyc1));
        dlzc = zeros(numel(lzc1));
    end
    % filter the values of xc and fill R
    Rx = filter_vals(Rx, dlxc, thresh, l, lab, lxc1, lxc2);
    Ry = filter_vals(Ry, dlyc, thresh, l, lab, lyc1, lyc2);
    Rz = filter_vals(Rz, dlzc, thresh, l, lab, lzc1, lzc2);
    %     if(maxxc>0) lxc1 = lxc1/maxxc;end
    %     if(maxyc>0) lyc1 = lyc1/maxyc;end
    %     if(maxzc>0) lzc1 = lzc1/maxzc;end
    %     if(maxxc>0) lxc2 = lxc2/maxxc;end
    %     if(maxyc>0) lyc2 = lyc2/maxyc;end
    %     if(maxzc>0) lzc2 = lzc2/maxzc;end
    
    if plotflag
        dfig(1);h = subplot(numel(0:L_max),1,plot_count);plot(dlxc,'--rs','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        set(h,'XTick',1:numel(lab));set(h,'XTickLabel',lab);
        dfig(2);h = subplot(numel(0:L_max),1,plot_count);plot(dlyc,'--rs','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        set(h,'XTick',1:numel(lab));set(h,'XTickLabel',lab);
        dfig(3);h = subplot(numel(0:L_max),1,plot_count);plot(dlzc,'--rs','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        set(h,'XTick',1:numel(lab));set(h,'XTickLabel',lab);
        plot_count = plot_count + 1;
    end
end
display_R(Rx);
display_R(Ry);
display_R(Rz);

%%
function display_R(R)
if numel(R)
f1 = figure('Position',[100 100 1000 900]);
cnames = {'Coefficient label','L','K', 'Score','index', 'CLK_o', 'CLK'};
t3 = uitable('Parent',f1,'Data',R,'ColumnName',cnames,...
    'Position',[20 20 960 860]);set(t3,'ColumnWidth',{120});
else
    disp('empty vector --- skipping');
end

function R = filter_vals(R,dlxc, thresh, l, lab, lxc1, lxc2)
cnt = 1;
for m = -l:l
    if dlxc(cnt)>thresh,
        R(end+1,1) = {lab{cnt}};
        R(end,2) = {l};
        R(end,3) = {m};
        R(end,4) = {dlxc(cnt)};
        R(end,5) = {(l-1 + 1)^2 + cnt};
        R(end,6) = {lxc1(cnt)};
        R(end,7) = {lxc2(cnt)};
    end
    cnt = cnt + 1;
end

function lab = label_gen(L)
lab= {};
count = 1;
for m = -L:L
    str = '';
    if (m<0),  str = [num2str(L) '_m' num2str(abs(m))];
    else
        str = [num2str(L) '_' num2str(m)];
    end
    lab{count} = str;
    count = count + 1;
end












