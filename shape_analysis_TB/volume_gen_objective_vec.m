function [Fun_vec] = volume_gen_objective_vec(Xo)
%%% For use with fgoalattain
%%% Objective function that calculates the residual vector for the difference between the 
%%% true raw data volume and the convolved  model volume. the vector Xo has the clks, then the shifts
global Y_LK PSF_fft x_pixel y_pixel z_pixel RAW F normXo plotflag fac frame_vec Iobj chisq ind_RAW
global ind_parm full_parm wbvec_x wbvec_y wbvec_z

Xo = Xo.*normXo;
clks = full_parm;
clks(ind_parm) = Xo(1:end-3);
clks(ind_parm) = Xo;
% xshift = Xo(end-2);yshift = Xo(end-1);zshift = Xo(end);
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc);zclks =clks(2*nc+1:end);
intensity_gen(xclks, yclks, zclks); % generate the 3D intensity image Iobj from the clks
Iobj = Iobj.*(Iobj(:)\RAW(:));    % % scale by using linear least squares
% chisq(:) = 0;
chisq = sum((RAW(ind_RAW)-Iobj(ind_RAW)).^2);   % The CHISQUARE functional
% [Bx By Bz] =  single_mode_bending_energy(xclks(:),yclks(:),zclks(:));                              % The smoothing functional
[B] =  bending_energy(xclks(:),yclks(:),zclks(:));                              % The smoothing functional

% Fun_vec = [chisq B Bx By Bz]';
Bx = wbvec_x(:)'*xclks(:);By = wbvec_y(:)'*yclks(:);Bz = wbvec_z(:)'*zclks(:);
Fun_vec = [chisq B abs(Bx) abs(By) abs(Bz)];

%for ix = 1:size(Iobj,4),imraw = RAW(:,:,1,ix);imobj = Iobj(:,:,1,ix);R = [R; (imraw(:)-imobj(:))];end

global itercount;
if plotflag,
    itercount = itercount + 1;%disp(itercount);disp( sum(R.^2)/(length(R)-length(Xo)));
    if mod(itercount,10)==0
    str = sprintf('Iteration: %d     Chisq: %.4f    Bx: %.4f By: %.4f Bz: %.4f',itercount, Fun_vec(1), Fun_vec(2), Fun_vec(3), Fun_vec(4));
    cla;montage(mat2gray(Iobj));title(str);drawnow;
    end
end