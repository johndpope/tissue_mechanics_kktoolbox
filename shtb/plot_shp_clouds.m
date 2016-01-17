function plot_shp_clouds(SHP)
%%% Input: SHP cell array, with elements that are matrices T, where each
%%% row in T represents one shape. 


for ix  = 1:length(SHP),
 rand('seed',402);   
    dfig;
     plot_shps(SHP{ix}, 2, 'random');view(90,0);xlabel('x');ylabel('y');zlabel('z');
lighting phong;camlight;axis equal;
%     plot_shps(cs2nocs(SHP{ix}));
end

