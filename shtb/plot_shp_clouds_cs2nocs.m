function plot_shp_clouds_cs2nocs(SHP)
%%% Input: SHP cell array, with elements that are matrices T, where each
%%% row in T represents one shape. 
figure;
rand('seed',402);
for ix  = 1:length(SHP),
     plot_shps(cs2nocs(SHP{ix}));
end

view(3);xlabel('x');ylabel('y');
lighting phong;camlight
