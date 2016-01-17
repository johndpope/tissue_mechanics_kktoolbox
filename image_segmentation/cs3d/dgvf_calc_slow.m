function [un, vn, wn, Ix, Iy, Iz, x, y, z] = dgvf_calc(I, niter, miu, dt, dx, dy, dz, distributed)
%% generate the diffusion gradient vector field as in Xu and Prince 1998
verbose = 0;
[Ix, Iy, Iz] = gradient(double(I));
[x,y, z] = meshgrid(1:size(I,2), 1:size(I,1), 1:size(I,3));
r = miu*dt/dx/dy;
F = sqrt(Ix.^2 + Iy.^2 + Iz.^2); % Energy E2 as in Xu and Prince 1998
if (dt>(dx*dy*dz/6/miu)), disp('convergence not guaranteed!!!! Resetting dt');dt =(dx*dy*dz/6/miu/2);disp(dt);end
[fx, fy, fz] = gradient(F);
f = sqrt(fx.^2 + fy.^2 + fz.^2);
b = fx.^2 + fy.^2 + fz.^2;
c1 = b.*fx;
c2 = b.*fy;
c3 = b.*fz;

disp('Calculating 3D diffusion gradient vector field');
if verbose,tic;end
if distributed,
    %%%%%%%%%%% use three processors to do this %%%%%%%%%%%%%%%%
    inp1{1} = Ix;inp1{2} = niter;inp1{3} = b;inp1{4} = r;inp1{5} = dt;inp1{6} = c1;
    inp2{1} = Iy;inp2{2} = niter;inp2{3} = b;inp2{4} = r;inp2{5} = dt;inp2{6} = c2;
    inp3{1} = Iy;inp3{2} = niter;inp3{3} = b;inp3{4} = r;inp3{5} = dt;inp3{6} = c3;
    sched = findResource('scheduler', 'type', 'local');
    set(sched, 'DataLocation','./')
    job1 = createJob(sched);
    set(job1,'FileDependencies',{'ds.m'})
    warning off;createTask(job1, @ds, 1, { {Ix,niter,b,r,dt,c1} {Iy,niter,b,r,dt,c2} {Iz,niter,b,r,dt,c3}});warning on;
    submit(job1);
    waitForState(job1);
    results = getAllOutputArguments(job1)';disp('done!');
    errmsgs = get(job1.Tasks, {'ErrorMessage'});
    nonempty = ~cellfun(@isempty, errmsgs);
    celldisp(errmsgs(nonempty));
    un = results{1};vn = results{2}; wn = results{3};
    %%%%%%%%% end of using the distributed computing %%%%%%%%%%%%%
else
    %%%%%%%%% alternatively, use old (non-distributed) method %%%%
    un = Ix;
    vn = Iy;
    wn = Iz;
    %% prepare initial quantities to avoid memory allocations later on
    %% and excessive calls to circshift
    unp1 = zeros(size(I));
    vnp1 = zeros(size(I));
    wnp1 = zeros(size(I));
    
    indx_mx = reshape(find(ones(size(I))), size(I));  % just get a matrix with the linear indices
    xm1 = circshift(indx_mx,[-1   0   0]);  % generate linear indices with the circular shifts
    x1  = circshift(indx_mx,[ 1   0   0]);
    ym1 = circshift(indx_mx,[ 0  -1   0]);
    y1  = circshift(indx_mx,[ 0   1   0]);
    zm1 = circshift(indx_mx,[ 0   0  -1]);
    z1  = circshift(indx_mx,[ 0   0   1]);

    
    for n = 1:niter,        %%loop over the desired time-steps (Xu Prince 1998 p.363
        if verbose,disp(['Iteration : ' num2str(n)]);end

        unp1(:) = (1-b(:).*dt).*un(:)+ ...
            r.*(  ...
            un(xm1(:)) + un(x1(:)) + un(ym1(:)) +...
            un(y1(:)) + un(zm1(:)) + un(z1(:))- 6.*un(:))...
            + c1(:).*dt;

        vnp1(:) = (1-b(:).*dt).*vn(:)+ ...
            r.*(  ...
            vn(xm1(:)) + vn(x1(:)) + vn(ym1(:)) +...
            vn(y1(:)) + vn(zm1(:)) + vn(z1(:))- 6.*vn(:))...
            + c1(:).*dt;

        wnp1(:) = (1-b(:).*dt).*wn(:)+ ...
            r.*(  ...
            wn(xm1(:)) + wn(x1(:)) + wn(ym1(:)) +...
            wn(y1(:)) + wn(zm1(:)) + wn(z1(:))- 6.*wn(:))...
            + c1(:).*dt;
%           
%         vnp1 = (1-b.*dt).*vn+ r.*(  ...
%               circshift(vn,[-1  0  0])...
%             + circshift(vn,[ 1  0  0])...
%             + circshift(vn,[ 0 -1  0] )...
%             + circshift(vn,[ 0  1  0] )...
%             + circshift(vn,[ 0  0 -1] )...
%             + circshift(vn,[ 0  0  1])...
%             - 6.*vn)...
%             + c2.*dt;
%         wnp1 = (1-b.*dt).*wn+ r.*(  ...
%               circshift(wn,[-1  0  0])...
%             + circshift(wn,[ 1  0  0])...
%             + circshift(wn,[ 0 -1  0] )...
%             + circshift(wn,[ 0  1  0] )...
%             + circshift(wn,[ 0  0 -1] )...
%             + circshift(wn,[ 0  0  1])...
%             - 6.*wn)...
%             + c3.*dt;
        un(:) = unp1(:);
        vn(:) = vnp1(:);
        wn(:) = wnp1(:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose,
toc
step = 6;
quiver3(x(1:step:end,1:step:end,1:step:end),...
        y(1:step:end,1:step:end,1:step:end),...
        z(1:step:end,1:step:end,1:step:end),...
        un(1:step:end,1:step:end,1:step:end),...
        vn(1:step:end,1:step:end,1:step:end),...
        wn(1:step:end,1:step:end,1:step:end));
drawnow;axis equal;view(3)
end