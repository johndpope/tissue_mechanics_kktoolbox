function [results3, com] = filter_objects(results2)
%%% filter rendered results
afac = 8;
vfac = 8;
vredfac = 8;

P.pre.vmin = 5;
P.pre.vmax = inf;
P.pre.amin = 5;
P.pre.amax = inf;
P.pre.vred_min = 0.05;
results3 = {};
for ix = 1:length(results2)
    fprintf('\nFiltering time point: %d\n', ix);
    all_X = results2{ix}.all_X;
    all_F = results2{ix}.all_F;
    all_X_out = results2{ix}.all_X_out;
    all_F_out = results2{ix}.all_F_out;
    %%% remove empty or invalid entries in all_X and all_X_out
    remix = [];
    for gix = 1:length(all_X_out)
        if isempty(all_X_out{gix}),remix = [remix gix];end
    end
    all_X_out(remix) = [];
    all_F_out(remix) = [];
    fprintf('Deleted %d invalid entries in all_X_out\n', length(remix));

    
    remix = [];
    for gix = 1:length(all_X)
        if isempty(all_X{gix}),remix = [remix gix];end
    end
    all_X(remix) = [];
    all_F(remix) = [];
    fprintf('Deleted %d invalid entries in all_X\n', length(remix));
    
    %% calculate the areas and volumes
    fprintf('\nCalculating geometrical properties of all objects\n');
    A = zeros(length(all_X_out),1);
    V = zeros(length(all_X_out),1);
    v = zeros(length(all_X_out),1);
    for gix = 1:length(all_X_out)
            [area vol red_vol] = triangulated_props(all_X_out{gix}, all_F_out{gix}, 0);
            A(gix) = area;
            V(gix) = vol;
            v(gix) = red_vol;
            %disp([area vol red_vol]);
    end
    %% filter the objects by area, volume, and reduced volume
    all_X2 = all_X_out;
    all_F2 = all_F_out;
    A2 = A;
    V2 = V;
    v2 = v;
    %%% remove those that don't lie within the limits
    indx = logical([V<P.pre.vmin | V>P.pre.vmax | v<P.pre.vred_min | A<P.pre.amin | A>P.pre.amax]);
    all_X2(indx) = [];
    all_F2(indx) = [];
    A2(indx) = [];
    V2(indx) = [];
    v2(indx) = [];
    disp('Enforcing constraints');
    fprintf(' ...found %d\n', length(indx));
    %%% remove outliers for area
    disp('Detecting and removing area outliers');
    indx = determine_outliers(A2, afac);
    fprintf(' ...found %d\n', length(indx));
    all_X2(indx) = [];
    all_F2(indx) = [];
    A2(indx) = [];
    V2(indx) = [];
    v2(indx) = [];
    %%% remove outliers for volume
    disp('Detecting and removing volume outliers');
    indx = determine_outliers(V2, vfac);
    fprintf(' ...found %d\n', length(indx));
    all_X2(indx) = [];
    all_F2(indx) = [];
    A2(indx) = [];
    V2(indx) = [];
    v2(indx) = [];
     %%% remove outliers for reduced volume
    disp('Detecting and removing reduced volume outliers');
    indx = determine_outliers(v2, vredfac);
    fprintf(' ...found %d\n', length(indx));
    all_X2(indx) = [];
    all_F2(indx) = [];
    A2(indx) = [];
    V2(indx) = [];
    v2(indx) = [];

    X = center_of_mass_all_X(all_X2);
    
    results3{ix}.all_X_out = all_X_out;
    results3{ix}.all_F_out = all_F_out;
    results3{ix}.all_X2 = all_X2;
    results3{ix}.all_F2 = all_F2;
    results3{ix}.com = X;
    results3{ix}.A = A;
    results3{ix}.V = V;
    results3{ix}.v = v;
    results3{ix}.A2 = A2;
    results3{ix}.V2 = V2;
    results3{ix}.v2 = v2;
end
%% construct the com cell array
com = cell(length(results3),1);
for ix = 1:length(results3)
    com{ix} = results3{ix}.com;
end
% % %% plot the result
% % axis(gca);
% % for ix = 1:length(results3)
% %     cla;
% %     rand('seed', 100);
% %     disp(ix);
% %     A = results3{ix}.A2;
% %     V = results3{ix}.V2;
% %     v = results3{ix}.v2;
% %     
% %     c = colormap(jet(max(ceil(V))));
% %     %%%%%%%%%%% plot the objects
% %     for oix = 1:size(results3{ix}.all_X2,2),   % loop over the shapes
% %         A = results3{ix}.A2(oix);
% %         V = results3{ix}.V2(oix);
% %         v = results3{ix}.v2(oix);
% %         map_parm = V;
% %         X = results3{ix}.all_X2{oix};
% %         F = results3{ix}.all_F2{oix};
% %         %a = rand(1,3);      % random RGB color specification
% %         patch('Vertices', X, 'Faces', F,'FaceColor', c(round(map_parm),:),'EdgeColor','none','FaceAlpha',1);
% %         %drawnow;
% %         hold on
% %     end
% %     %plot_all_X(results3{ix}.all_X2, results3{ix}.all_F2, 'r');
% %     view(0, -90);axis equal;camlight; axis off; axis tight;drawnow
% % end






















