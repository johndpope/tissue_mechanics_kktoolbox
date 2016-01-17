function averageX = pd_montecarlo(mc_opts, A_o, q, v_o, a_o_bar, X_o, Y_LK_tri,t, p,  flag,lambda_shear, filename)
%%% calculate the global minimum of the energy returned by objfun_ade
%%% within the tolerances of constraints on the total surface area, and reduced volume

averageX = []; counter = 0; rejected_count = 0;
E_vec = []; Xsat = [];

mxfun = mc_opts.mxfun;                              % the number of shapes to be tested for constraint satisfaction at a time
target_sample_size = mc_opts.target_sample_size ;   % calculate batches of trials until this number of shapes satisfying the constraints is reached
trial_limit = mc_opts.trial_limit;                  % maximum number of batches to try
beta = mc_opts.beta;
lambda = mc_opts.lambda;
area_tol = mc_opts.area_tol;
vol_tol = mc_opts.vol_tol;


maxStore = target_sample_size * 10;
str = sprintf('Batch size = %d\n Target sample # = %d\n Max. trials = %d',...
    mxfun, target_sample_size, trial_limit); disp(str);

if mc_opts.precalc,
    warning off;precalc_constraints_vectorized;disp('finished memory allocation');
end
    while (size(Xsat,1)<target_sample_size && counter < trial_limit)
        [Currentsat, A, V] = satisfy_constraints_vec(X_o,A_o,Y_LK_tri,v_o,...
            lambda, mxfun, area_tol, vol_tol, mc_opts.miu);
        Xsat = [Xsat; Currentsat];
        counter = counter + 1;
        nshapes = size(Xsat,1);
%         str = sprintf('%d  shapes of %d  ok. Area range: %.2f to %.2f , V: %.2f to %.2f.',...
%             nshapes, counter * mxfun, min(A), max(A), min(V), max(V));disp(str)
                str = sprintf('%d  shapes of %d  ok. V: %.2f to %.2f.',...
            nshapes, counter * mxfun, min(V), max(V));disp(str)
    end
    save data_MC_shapes Xsat;


load data_MC_shapes;
counter = 0;rejected_count = 0;E_vec = [];
a_o_bar = a_o_bar(1);
v_o     = v_o(1);
A_o     = A_o(1);
nshapes = size(Xsat,1);
storeX   = zeros(nshapes,length(X_o));

if nshapes,
    if satisfy_constraints(X_o,A_o,Y_LK_tri,v_o,area_tol, vol_tol),
        [Eold] = objfun_ade(X_o,A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear,flag);
        Xaccepted = X_o;
    else
        disp('\n Old shape does not satisfy constraints...picking from new set./n');
        [Eold] = objfun_ade(Xsat(1,:),A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear,flag);
        Xaccepted = Xsat(1,:);
    end
    counter = counter + 1;
    storeX(counter,:) = Xaccepted;
    for ix  = 1:nshapes,
        if mod(ix,1000)==0,disp(ix);end
        Xtry = Xsat(ix,:);
        [Enew] = objfun_ade(Xtry,A_o, q,Y_LK_tri, v_o,a_o_bar, lambda_shear,flag);
        if Enew<Eold,               % accept 
            counter = counter + 1;
            E_vec = [E_vec Enew];
            disp([ix 0]);
            storeX(counter,:) = Xtry;
            Xaccepted = Xtry;
            Eold = Enew;
            if counter == maxStore, break;end;
        elseif (exp(-beta*(Enew-Eold)))>rand,   % accept (Metropolis +ve)
            counter = counter + 1;
            E_vec = [E_vec Enew];
            disp([ix 1]);
            storeX(counter,:) = Xtry;
            Xaccepted = Xtry;
            Eold = Enew;
            if counter == maxStore, break;end
%         else,                                   % reject
%             counter = counter + 1;
%             rejected_count = rejected_count + 1;
%             storeX(counter,:) = Xaccepted;
%             if counter == maxStore, break;end
        end
    end
    if counter,
%        averageX = Xaccepted;
        averageX = sum(storeX)./counter;
        [E,v_red,A, V] = objfun_ade(averageX,A_o, q,Y_LK_tri, v_o,a_o_bar,lambda_shear, 2);axis equal;
%         figure;plot(E_vec);title('Progress of the Shape energy during the MC simulation');
%         disp([num2str(rejected_count) ' rejected']);
    end
else averageX = [];
    
end


%%%%%%%%%%%%% select based on Area and volume
function [ok] = satisfy_constraints(clks,A_o,Y_LK,v_o, Areatol, Voltol)
ok = 0;
nc = round(length(clks)/3);xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
X = [Y_LK*xclks(:)  Y_LK*yclks(:) Y_LK*zclks(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculates the Area(A), the volume (V) and the local mean and total mean curvatures (H and h) of the shape.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global u v w crossqpr C
%% calculation of the area and volume
u(:) = X(:,1); v(:) = X(:,2); w(:) = X(:,3);
crossqpr(:) = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],...
    [u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));   %
A = sum(twoA)/2; % this is the total area
n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));
Vo = 4/3*pi*(A/4/pi)^(3/2);v_red = V/Vo;
v_red_o = v_o./(4/3*pi*(A_o/4/pi).^(3/2));
if abs(v_red-v_red_o)<=Voltol,ok = 1;end

    
    