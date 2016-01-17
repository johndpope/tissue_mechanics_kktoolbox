function [R] = volume_gen_objective(Xo)
%%% Objective function that calculates the residual vector for the difference between the 
%%% true raw data volume and the convolved  model volume. the vector Xo has the clks, then the shifts
global Y_LK PSF_fft x_pixel y_pixel z_pixel RAW F normXo plotflag fac frame_vec Iobj R ind_RAW
global ind_parm full_parm alpha

Xo = Xo.*normXo;
clks = full_parm;
clks(ind_parm) = Xo(1:end-3);
%obj_scale = Xo(end-3);x_pixel = x_pixel *obj_scale; y_pixel = y_pixel *obj_scale;z_pixel = z_pixel * obj_scale;

xshift = Xo(end-2);yshift = Xo(end-1);zshift = Xo(end);

nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
intensity_gen(xclks, yclks, zclks, xshift, yshift, zshift); % generate the 3D intensity image Iobj from the clks
% scale by using linear least squares (should be done in one step)
%%Iobj = Iobj.*(Iobj(:)\RAW(:));    % determine scale using the QR decomposition
for ix = 1:size(RAW,4),    %loop over the frames
    imIobj = Iobj(:,:,1,ix);imRAW = RAW(:,:,1,ix);
    A = imIobj(:);b = imRAW(:); scale = A\b;    % determine x using the QR decomposition
    s(ix) = scale;
    Iobj(:,:,1,ix) = Iobj(:,:,1,ix).*scale;
end
%%% Calculate the residual
R(:) = 0;
% for ix = 1:size(Iobj,4), 
%     imraw = imresize(RAW(:,:,1,ix), fac);
%     imobj = imresize(Iobj(:,:,1,ix), fac);
%     R = [R; (imraw(:)-imobj(:))];
% end
% imobj = Iobj(:,:,1,frame_vec);R(:) = (data_vec-imobj(:));
R = sum((RAW(ind_RAW)-Iobj(ind_RAW)).^2);
reg = Xo.^-2;reg(isinf(reg))=999999999;
REG =  sum(reg);
R = [R REG]';

%for ix = 1:size(Iobj,4),imraw = RAW(:,:,1,ix);imobj = Iobj(:,:,1,ix);R = [R; (imraw(:)-imobj(:))];end

global itercount;
if plotflag,
    itercount = itercount + 1;%disp(itercount);disp( sum(R.^2)/(length(R)-length(Xo)));
    if mod(itercount,10)==0
    str = sprintf('Iteration: %d     R: %.4f    reg: %.4f',itercount, R(1), R(2));
    cla;montage(mat2gray(Iobj));title(str);drawnow;
    end
end