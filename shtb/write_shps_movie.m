function M = write_shps_movie(S, gdim, nframes, filename)
%%% takes a series of SHP shapes and rotates them in a figure recording the intermediate frames
%%% Output is an mpeg movie and two sets of movie frames in subdirectories ./tifs and ./jpgs
%%% Input: S is a nxm matrix where n is the number of SHP objects and m is the lenght of coefficients
%%% vector (just like X_o).
%%% Example: M = write_shps_movie(S, 60, 100, 'cyst03');
%%% Copyright: Khaled Khairy, January 2009
%%% Also see: write_shp_movie, plot_shp_shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, gdim = 120; nframes = 30;end
if nargin<3, nframes = 30;end
if nargin<4, filename = 'untitled';end
figure(1);
g = 0;dg = deg2rad(360/nframes);
X = [];     % this is where we will store all coordinates just in one matrix
Xsmx = [];  % store all coordinates individually for each shape
Csmx = [];  % store all connectivities for each shape

if mod(size(S,2),3)==1, S = S(:,1:end-1);end
Lmax = sqrt(size(S,2)/3)-1;
cla;
%%% determine center of mass
disp('Calculating center of mass...');
for six = 1:size(S,1),
    [Xs,Cs]=sh_calc(S(six,:), gdim, Lmax);
    X = [X ; Xs];
    Xsmx(:,:,six) = Xs;
    Csmx(:,:,six) = Cs;
end
cm = sum(X,1)/length(X);
disp('Done!');


%%% plot all shps
for six = 1:size(S,1),
    [Xs,Cs]=sh_calc(S(six,:), gdim, Lmax);
    Xs(:,1) = Xs(:,1)-cm(1);          % to center around the center of mass
    Xs(:,2) = Xs(:,2)-cm(2);          % to center around the center of mass
    Xs(:,3) = Xs(:,3)-cm(3);          % to center around the center of mass
    Xsmx(:,:,six) = Xs;               % store all coordinates after translation
    patch('Vertices', Xs, 'Faces', Cs,'FaceColor', 'r', 'EdgeColor', 'none');
end
view([-180 -84]);
daspect([1 1 0.25]);
axis tight;
xlim([-max(X(:,1))*fac max(X(:,1))*fac]);
ylim([-max(X(:,2))*fac max(X(:,2))*fac]);
zlim([-max(X(:,3))*fac max(X(:,3))*fac]);
drawnow;

M = moviein(nframes);
X_start = X;
if exist([pwd '/tifs'])~= 7, mkdir([pwd '/tifs']);end
if exist([pwd '/jpgs'])~= 7, mkdir([pwd '/jpgs']);end

save temp
for ix = 1:nframes,
    
    disp(['Generating frame ' num2str(ix)]);
    cla;
    %%% plot all shps
    for six = 1:size(S,1),
        Xss = rotate_cartesian_vectors(Xsmx(:,:,six),[0 g 0]);
        
        
        % if you want several colors mapped to volume uncomment below
        % %         u = Xss(:,1); v = Xss(:,2); w = Xss(:,3);
        % %         crossqpr = cross([u(Csmx(:,2))-u(Csmx(:,1)) v(Csmx(:,2))-v(Csmx(:,1)) w(Csmx(:,2))-w(Csmx(:,1))],...
        % %             [u(Csmx(:,3))-u(Csmx(:,1)) v(Csmx(:,3))-v(Csmx(:,1)) w(Csmx(:,3))-w(Csmx(:,1))],2);
        % %         twoA = sqrt(sum(crossqpr.*crossqpr,2));Area = sum(twoA)/2;
        % %         F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
        % %         V = -sum(1/3*(dot(n,[u(Csmx(:,1)) v(Csmx(:,1)) w(Csmx(:,1))], 2).*twoA./2));
        % %         Vo = 4/3*pi*(Area/4/pi)^(3/2);v_red = V/Vo;
        % %         V = abs(V);v_red = abs(v_red);
        a = rand(1,3);      % random RGB color specification
        if sum(a)> 2.8; a = [1 0 1];end
        patch('Vertices', Xss, 'Faces', Csmx(:,:,six),'FaceColor', a,'EdgeColor','none','FaceAlpha', 1);
        %%%  if you want one color uncomment below
        %         patch('Vertices', Xss, 'Faces', Csmx(:,:,six),'FaceColor', 'r', 'EdgeColor', 'none');
    end
    light('Position',[0.9446 -0.4607 -0.9463]);
    %light('Position',[-1.393 0.2456 0]);
    light('Position',[-0.2881 -0.3613 0.8868]);
    %view([33 24]);
    %view([0 -90]);
    view([23 26]);
    daspect([ 1 1 0.25]);
    %lighting gouraud;
    lighting phong;
    axis off;
    drawnow;
    
    %%%%%% store the frame
    directory = [pwd '/tifs'];
    if ix<10
        number = sprintf('00%d',ix);
    elseif ix<100
        number = sprintf('0%d',ix);
    elseif ix<1000
        number = sprintf('%d',ix);
    end
    
    M(:,ix) = getframe;
    Image = getframe;
    P = frame2im(Image);    % Convert to a image representation that  Matlab can handle
    extension = 'tif';
    framename = sprintf('%s/%s.%s', directory, [filename number],extension);
    %imwrite(P,framename, extension) % Finally write individual images
    str = sprintf('print -dtiff -r300 %s', framename);eval(str);
    
    
    directory = [pwd '/jpgs'];
    extension = 'jpg';     % for a directory called "images"
    framename = sprintf('%s/%s.%s', directory, [filename number],extension);
    %     print ('-dppmraw', framename);
    str = sprintf('print -djpeg -r300 %s', framename);eval(str);
    
    framename = sprintf('%s/%s.%s', directory, [filename number],'ppm');
    str = sprintf('print -dppmraw -r300 %s', framename);eval(str);
    %%%%%% update the rotation
    g = g + dg;
    
    %zoom(1/zoomfac);
end
zoom(zoomfac);drawnow
%crop_frames(directory, filename, extension, nframes, [300 900], [150 650]);

save movie M;


% % %%% encode the movie to avi if on linux and with mencoder installed
% % cd jpgs;
% % str = sprintf('!mencoder "mf://*.jpg" -mf fps=10 -o %s_movie.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800;',...
% %     filename);eval(str);
% % 
% % str = sprintf('!ffmpeg -r 10 -y -b 2000 -i %s%%03d.jpg %s_movie_from_jpg.mp4', filename, filename);eval(str);
% % str = sprintf('!ffmpeg -r 10 -y -b 2000 -i %s%%03d.ppm %s_movie_from_ppm.mp4', filename, filename);eval(str);
% % cd ..


end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_rot = rotate_cartesian_vectors(X,ang)
%% rotate the shape  X by the Euler angles a b g (radians).
a = ang(1);b = ang(2);g = ang(3);
%% rotation conventions are y-z-y
Rg = rot_mx(g,2);
Rb = rot_mx(b,1);
Ra = rot_mx(a,2);
X_rot = [Ra*Rb*Rg*X']';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,C, Y_LK]=sh_calc(X_o, gdim, Lmax, X, C, Y_LK)
X_o = tr(X_o, Lmax);
nc = length(X_o)/3;
xclks = X_o(1:nc);
yclks = X_o(nc+1:2*nc);
zclks = X_o(2*nc+1:end);
if nargin<4,
    P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
    [t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
    X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
    Y_LK = get_basis(t',p',gdim,Lmax);
end
X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crop_frames(directory, prfx, extension, nframes, x, y)
here = pwd;
cd(directory);
for ix = 1:nframes,
    if ix<10
        str = sprintf('%s00%d.%s',prfx, ix,extension);
    elseif ix<100
        str = sprintf('%s0%d.%s',prfx,ix, extension);
    elseif ix<1000
        str = sprintf('%s%d.%s',prfx, ix, extension);
    end
    im = imread(str, extension);
    imcropped = im(y(1):y(2), x(1):x(2));
    imwrite(imcropped, str, extension);
end
cd ..
end







