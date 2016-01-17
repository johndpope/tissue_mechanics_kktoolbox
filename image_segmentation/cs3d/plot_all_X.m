function plot_all_X(all_X, all_F,c)
% plot the objects represented by the cell arrays all_X and all_F
%     clf;axis on;axis equal;
%     view(3);
%     light;
%     lighting gouraud;
hold on;
%daspect([1 1 1]);
if nargin ==2,
    for ix = 1:size(all_X,2),   % loop over the shapes
        X = all_X{ix};
        F = all_F{ix};
        a = rand(1,3);      % random RGB color specification
        patch('Vertices', X, 'Faces', F,'FaceColor', (a),'EdgeColor','none','FaceAlpha',1);
        %drawnow;
        hold on
    end
    %hold off
else
        for ix = 1:size(all_X,2),   % loop over the shapes
        X = all_X{ix};
        F = all_F{ix};
        %a = rand(1,3);      % random RGB color specification
        patch('Vertices', X, 'Faces', F,'FaceColor', c,'EdgeColor','none','FaceAlpha',1);
        %drawnow;
        hold on
    end
    %hold off
end