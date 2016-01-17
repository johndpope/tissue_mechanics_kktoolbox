function M = write_shp_movie(X_o, gdim, nframes, filename)
%%% takes an SHP shape and rotates it in a figure recording the intermediate frames
%%% Output is an mpeg movie and two sets of movie frames in subdirectories ./tifs and ./jpgs
% for example: M = write_shp_movie(X_o, 60, 100, 'Mushroom');
% Author: Khaled Khairy, January 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, gdim = 120, nframes = 30;end
if nargin<3, nframes = 30;end
if nargin<4, filename = 'untitled';end
figure;
[xc, yc, zc] = get_xyz_clks(X_o);
xc(1) = 0;yc(1) = 0;zc(1) = 0;
X_o = [xc(:)' yc(:)' zc(:)'];   % just to center at zero

g = 0;dg = deg2rad(360/nframes);

Lmax = sqrt(length(X_o)/3)-1;
[X,C]=sh_calc(X_o, gdim, Lmax);
cm = sum(X,1)/length(X);
X(:,1) = X(:,1)-cm(1);          % to center around the center of mass
X(:,2) = X(:,2)-cm(2);          % to center around the center of mass
X(:,3) = X(:,3)-cm(3);          % to center around the center of mass

cla;patch('Vertices', X, 'Faces', C,'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
axis tight;
fac = 5;
xlim([-max(X(:,1))*fac max(X(:,1))*fac]);
ylim([-max(X(:,2))*fac max(X(:,2))*fac]);
zlim([-max(X(:,3))*fac max(X(:,3))*fac]);

M = moviein(nframes);
X_start = X;
if exist([pwd '/tifs'])~= 7, mkdir([pwd '/tifs']);end
if exist([pwd '/jpgs'])~= 7, mkdir([pwd '/jpgs']);end
for ix = 1:nframes,

    cla;patch('Vertices', X, 'Faces', C,'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    light('Position',[0.9446 -0.4607 -0.9463]);
    light('Position',[-1.393 0.2456 0]);
    light('Position',[-0.2881 -0.3613 0.8868]);
    view([33 24]);
    daspect([1 1 1]);
    lighting gouraud;
    axis off;
    zoom(4);
    drawnow;

    %%%%%% store the frame
    M(:,ix) = getframe;
    Image = getframe;
    P = frame2im(Image);    % Convert to a image representation that  Matlab can handle
    directory = [pwd '/tifs'];

    if ix<10
        number = sprintf('00%d',ix);
    elseif ix<100
        number = sprintf('0%d',ix);
    elseif ix<1000
        number = sprintf('%d',ix);
    end


    extension = 'tif';write_shp_movie(X_o,40, 100);
    framename = sprintf('%s/%s.%s', directory, [filename number],extension);
    imwrite(P,framename, extension) % Finally write individual images

    directory = [pwd '/jpgs'];
    extension = 'jpg';     % for a directory called "images"
    framename = sprintf('%s/%s.%s', directory, [filename number],extension);
    %     print ('-dppmraw', framename);
    print ('-djpeg', framename);
    %%%%%% update the rotation
    g = g + dg;
    [X]=rotate_cartesian_vectors(X_start,[-g g g]);
    zoom(0.25);
end

crop_frames(directory, filename, extension, nframes, [300 900], [150 650]);

save movie M;
%%% encode the movie to avi if on linux and with mencoder installed
cd jpgs;
str = sprintf('!mencoder "mf://*.jpg" -mf fps=10 -o %S_movie.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800;',...
    filename);eval(str);

%str = sprintf('!ffmpeg -r 10 -b 1800 -i %s%%03d.jpg %s_movie.mp4', filename, filename);eval(str);
cd ..
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







