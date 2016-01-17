function [mosaic num] = gft_calc(un, vn, wn, x, y, z, thresh)
%%%% Gradient flow tracking 
tic

%%%%%%%%%%%%%% vectorized version %%%%%%%%%%%%%%%%%%%
disp('Gradient flow tracking (vectorized version)!');
xs = x;ys = y;zs = z;
xsp = x;ysp = y;zsp = z;
maxxr = max(x(:));maxyr = max(y(:));maxzr = max(z(:));
Ix = sub2ind(size(x), y(:), x(:), z(:));
counter = 0;
while counter<100,
   counter = counter + 1;disp(counter);
   Vx = sub2ind(size(x), ys(:), xs(:),zs(:));
   NV = sqrt((un(Vx).^2+vn(Vx).^2+wn(Vx).^2));
%    disp([x(vv) y(vv) z(vv) ;...
%        xs(vv) ys(vv) zs(vv);...
%        xsp(vv) ysp(vv) zsp(vv);...
%        un(Vx(vv)) vn(Vx(vv)) wn(Vx(vv))]);
   %un(20,40,35) vn(20,40,35) wn(20,40,35)]);
   xsp(Ix) = xs(Ix)+ round(un(Vx)./NV);xsp(xsp>=maxxr) = xs(xsp>=maxxr);xsp(xsp<=0) = 1;
   ysp(Ix) = ys(Ix)+ round(vn(Vx)./NV);ysp(ysp>=maxyr) = ys(ysp>=maxyr);ysp(ysp<=0) = 1;
   zsp(Ix) = zs(Ix)+ round(wn(Vx)./NV);zsp(zsp>=maxzr) = zs(zsp>=maxzr);zsp(zsp<=0) = 1;
   xs(:) = xsp(:);
   ys(:) = ysp(:);
   zs(:) = zsp(:);
   %hold on;plot3(xs(vv), ys(vv), zs(vv), '*k');drawnow;hold on;
end
step = 30;
hold on;plot3(xs(1:step:end), ys(1:step:end), zs(1:step:end), '*k');
% % x = x(:);y = y(:);z = z(:);
% % un = un(:); vn = vn(:); wn = wn(:);
% % sinkx = x(:);sinky = y(:);sinkz = z(:);
% % maxxr = max(x(:));maxyr = max(y(:));maxzr = max(z(:));
% % %theta_check = zeros(size(sinkx), 'double');
% % cum_check   = ones(size(sinkx), 'double');
% % Iix = sub2ind(size(I), y(:), x(:), z(:));    % general linear index into I-type matrices
% % %sinkxp = sinkx;sinkyp = sinky;sinkzp = sinkz; %
% % NV = sqrt((un(:).^2+vn(:).^2+wn(:).^2));
% % sinkxp = sinkx(:) + round(un(:)./NV);sinkxp(sinkxp>=maxxr) = sinkx(sinkxp>=maxxr);sinkxp(sinkxp<=0) = 1;
% % sinkyp = sinky(:) + round(vn(:)./NV);sinkyp(sinkyp>=maxyr) = sinky(sinkyp>=maxyr);sinkyp(sinkyp<=0) = 1;
% % sinkzp = sinkz(:) + round(wn(:)./NV);sinkzp(sinkzp>=maxzr) = sinkz(sinkzp>=maxzr);sinkzp(sinkzp<=0) = 1;
% % theta_check = 1;
% % sink = sum(sinkx.^2 + sinky.^2 + sinkz.^2);
% % sinkp = sum(sinkxp.^2 + sinkyp.^2 + sinkzp.^2);
% % counter = 0;
% % vv = 272118;
% % while counter<100%(any(sink-sinkp)), % as long as voxels are still flowing
% %     counter = counter + 1; disp(counter);
% %     disp([x(vv) y(vv) z(vv) ;sinkx(vv) sinky(vv) sinkz(vv);sinkxp(vv) sinkyp(vv) sinkzp(vv); un(vv) vn(vv) wn(vv)]);
% %  
% %     NV = sqrt(un(sinkx).^2 + vn(sinky).^2 + wn(sinkz).^2);
% %  %   t1   = [un(sinkx) vn(sinky) wn(sinkz)]./[NV NV NV];
% %     NVp = sqrt(un(sinkxp).^2 + vn(sinkyp).^2 + wn(sinkzp).^2);
% %  %   t2   = [un(sinkxp) vn(sinkyp) wn(sinkzp)]./[NVp NVp NVp];
% %  %   theta_check = (acos(sum(t1.*t2,2))<pi/2);
% %  %   cum_check(theta_check==0) = 0;
% %     
% %     sinkx(:)  = sinkxp(:);
% %     sinky(:)  = sinkyp(:);
% %     sinkz(:)  = sinkzp(:);
% %     sinkxp(:) = sinkx + cum_check.*theta_check.*round(un(sinkx)./NVp);
% %     sinkyp(:) = sinky + cum_check.*theta_check.*round(vn(sinky)./NVp);
% %     sinkzp(:) = sinkz + cum_check.*theta_check.*round(wn(sinkz)./NVp);
% %     
% %     sink = sum(sinkx.^2 + sinky.^2 + sinkz.^2);
% %     sinkp = sum(sinkxp.^2 + sinkyp.^2 + sinkzp.^2);
% %     
% % 
% % disp([x(vv) y(vv) z(vv) ;sinkx(vv) sinky(vv) sinkz(vv);sinkxp(vv) sinkyp(vv) sinkzp(vv); un(vv) vn(vv) wn(vv)]);
% % hold on;plot3(sinkx(vv), sinky(vv), sinkz(vv), '*k');drawnow;hold on;
% % end
toc
disp('Done!');
warning off;

%%%%%%%tic
% % % % % %%%%%%%%%%%%% non-vectorized version
% % % % % V = mat2gray(sqrt(un.^2 + vn.^2 + wn.^2));
% % % % % index = find(V>thresh);
% % % % % xr = x(index);xr = xr(:);
% % % % % yr = y(index);yr = yr(:);
% % % % % zr = z(index);zr = zr(:);
% % % % % tobetrkd = zeros([size(x,2) size(x,1) size(x,3)], 'uint16');tobetrkd(index) = 1;
% % % % % sinkx = zeros(size(x),'uint16');
% % % % % sinky = zeros(size(x),'uint16');
% % % % % sinkz = zeros(size(x),'uint16');
% % % % % disp('Calculating gradient flow tracking.');
% % % % % disp(['Number of points to track: ' num2str(length(yr))]);drawnow;
% % % % % for ix = 1:length(xr)-1,      %loop over relevant points and track them
% % % % %         if mod(ix,5000)==0,disp([num2str(ix/length(yr)*100) '%']);end
% % % % %         counter = 0;
% % % % %         X = [xr(ix) yr(ix) zr(ix)];
% % % % %         if tobetrkd(X(1), X(2), X(3)), keep_going = 1; else keep_going = 0;end
% % % % %         
% % % % %         Xp = X + round([un(X(2),X(1),X(3)) vn(X(2),X(1),X(3)) wn(X(2),X(1),X(3))]...
% % % % %             /sqrt(sum([un(X(2), X(1),X(3)) vn(X(2),X(1),X(3)) wn(X(2),X(1),X(3))].^2)));
% % % % %         if ~(any(Xp<=0) || Xp(2)>=size(x,1) || Xp(1)>=size(x,2) || Xp(3)>=size(x,3))
% % % % %             onpath = [];
% % % % %             while (keep_going),  % track the point until angle changes by > 90 degrees
% % % % %                 counter = counter + 1;
% % % % %                 vec = [un(X(2),X(1),X(3)) vn(X(2),X(1),X(3)) wn(X(2),X(1),X(3))];
% % % % %                 vecp = [un(Xp(2),Xp(1),Xp(3)) vn(Xp(2),Xp(1),Xp(3)) wn(Xp(2),Xp(1),Xp(3))];
% % % % %                 t1 = vec/sqrt(sum(vec.^2));
% % % % %                 t2 = vecp/sqrt(sum(vecp.^2));
% % % % % 
% % % % %                 if (isnan(t1./t2)), theta=pi;
% % % % %                 else
% % % % %                     theta = acos(sum(t1.*t2));
% % % % %                     %theta =acos(dot(t1,t2));    % pfi
% % % % %                 end
% % % % %                 if theta>pi/2,
% % % % %                     keep_going =0;
% % % % %                     sinkx(yr(ix), xr(ix), zr(ix)) = X(1); 
% % % % %                     sinky(yr(ix), xr(ix), zr(ix)) = X(2);
% % % % %                     sinkz(yr(ix), xr(ix), zr(ix)) = X(3);
% % % % %                     if ~isempty(onpath),
% % % % %                         sinkx(onpath(:,2), onpath(:,1),onpath(:,3)) = X(1);
% % % % %                         sinky(onpath(:,2), onpath(:,1),onpath(:,3)) = X(2);
% % % % %                         sinkz(onpath(:,2), onpath(:,1),onpath(:,3)) = X(3);
% % % % %                     end
% % % % %                 else
% % % % %                     tobetrkd(X(1), X(2), X(3))=0;
% % % % %                     onpath = [onpath; [X(1) X(2) X(3)]];
% % % % %                     %disp(['counter: ' num2str(counter) '   ' num2str([X]) '   ' num2str(rad2deg(theta))]);
% % % % %                     X = Xp;
% % % % %                     Xp = X + round([un(X(2),X(1),X(3)) vn(X(2),X(1),X(3)) wn(X(2),X(1),X(3))]...
% % % % %                         /sqrt(sum([un(X(2),X(1),X(3)) vn(X(2),X(1),X(3)) wn(X(2),X(1),X(3))].^2)));
% % % % %                     if any(Xp<=0) || Xp(2)>=size(x,1) || Xp(1)>=size(x,2) || Xp(3)>=size(x,3), keep_going = 0;end
% % % % %                     %hold on;plot3(Xp(1), Xp(2), Xp(3), '*k');drawnow;
% % % % %                 end
% % % % %             end
% % % % %         end
% % % % % end
% % % % % toc
% % % % % warning on;
% % % % % allix = find((sinkx+sinky+sinkz)>0);    %allix ->linear indices of voxels that have been tracked
% % % % % trkend =[sinkx(allix) sinky(allix) sinkz(allix)]; %trkend(trkend(:,1)==0) = []; % trkend-> three vector of target sink coordinates for each followed voxel
% % % % % hold on;plot3(trkend(:,1), trkend(:,2), trkend(:,3), '*k');axis equal;
% % % % % 
% % % % % 
% % % % % %% let's analyse the results and determine the number of sinks 
% % % % % 
% % % % % disp('Identifying and fusing sinks.');
% % % % % mosaic = zeros(size(x));
% % % % % for ix = 1:length(trkend),
% % % % %     mosaic(trkend(ix,2), trkend(ix,1), trkend(ix,3)) = 1;
% % % % % end
% % % % % L = bwlabeln(mosaic, 26);
% % % % % num = max(L(:));
% % % % % disp([' Initial number of sinks found is: ' num2str(num)]);
% % % % % figure;kk_montage(mat2gray((mosaic)));
% % % % % 
% % % % % disp('Assigning sinks to voxels');
% % % % % for six = 1:num,        % loop over the sinks and record the relevant coordinates
% % % % %     disp(['Processing object number: ' num2str(six)]);
% % % % %     vec = find(L==six);
% % % % %     [xs ys zs] = ind2sub(size(L),vec);
% % % % %     for ix = 1:length(allix),
% % % % %         for ixs = 1:length(vec)
% % % % %             if sinkx(allix(ix)) == ys(ixs)&& sinky(allix(ix))==xs(ixs)...
% % % % %                     && sinkz(allix(ix))== zs(ixs),
% % % % %                 mosaic(allix(ix)) = six;
% % % % %             end
% % % % %         end
% % % % %     end
% % % % % end
% % % % % disp('done!');


