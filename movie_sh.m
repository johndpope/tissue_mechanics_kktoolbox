function m = movie_sh(CLKS, gdim, da_o, im_sequence)
% make a movie from the shape sequence contained in CLKS
% CLKS is a nframes x length(clks) matrix
m = [];
if nargin==1, gdim = 40;end
currdir = pwd;
if im_sequence, cd('./movie_frames');end
for ix = 1:size(CLKS,1),
    disp(ix);
    plot_sh(CLKS(ix,:),gdim);
    %%% set up appearance of everything
    str = sprintf('\nda_o: %.5f   A: %.1f \\mum^2  V: %.1f \\mum^3',da_o(ix), 143, 100);
    title(str, 'fontname','Times', 'fontunits','points', 'fontsize', 16, 'fontweight', 'normal');
%     textobj = findobj('type', 'text');set(textobj, 'fontunits', 'points');    set(textobj, 'fontsize', 16);set(textobj, 'fontweight', 'normal');    set(textobj, 'fontname', 'Times');
    xlim([-5 5]);ylim([-5 5]);zlim([-2 2]);
%    view(3);
    view(115,35);
     lighting phong; drawnow;
%     m(ix) = getframe(gcf);
    
    if im_sequence,     % write image files for the movie
        
        if     ix <10,  str = sprintf('0_image_000%d.tif', ix);
        elseif ix <100, str = sprintf('0_image_00%d.tif',ix);
        elseif ix <1000,str = sprintf('0_image_0%d.tif',ix);
        end
        printstr = sprintf('print -dtiff -zbuffer -r600 %s;',str);
        eval(printstr);
    end
end
cd(currdir);
% disp('Writing movie to disc');
% movie2avi(m,'moviefile.avi', 'compression','Cinepak','quality',100, 'fps', 30);
% movie(m,1);