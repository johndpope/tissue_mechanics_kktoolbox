function [t,p, shear, area_res] = optimize_set(F, S,t, p,F_areas_relative,maxiter, DMflag,JacPatCalc, AreaCalcEqual)

global pcount 
pcount = 1;
flag = 1;
p = mod(p,2*pi);
t1 = t(1);p1 = p(1);%t2 = t(end); p2 = p(end);
%t = t(2:end-1);p = p(2:end-1);Xfixed = [t1 p1;t2 p2];
%t = t(2:end);p = p(2:end);
Xo = [t;p];nvert = length(t);
DM = []; invDM = [];F_areas_o = [];
if DMflag, [DM, invDM, F_areas_o] = triangle_deformation_init(F, S(:,1), S(:,2), S(:,3));end
% % % %% determine sparsity pattern for jacobian in order to be able to use the large-scale algorithm
if JacPatCalc
    warning off;disp('Calculating Jacobian Pattern'); nvert = length(Xo)/2;
    JacPat = sparse(length(F_areas_relative),nvert);
    for ix = 1:length(F),
        verts = F(ix,:);
        for vert = 1:length(verts),
            JacPat(ix,verts(vert)) = 1;
        end
    end
    %JacPat(:,1) = [];
    JacPat = sparse([JacPat JacPat]);
    JacPat = sparse([JacPat;zeros(4,size(JacPat,2))]);
    save JacPat JacPat
end

disp('Beginning optimization');
load JacPat;
%%%%%%%%%% Prepare quantities
global u v w crossqpr twoA F_areas
u = zeros(length(Xo)/2,1,'double');
v = zeros(length(Xo)/2,1,'double');
w = zeros(length(Xo)/2,1,'double');
crossqpr = zeros(size(F),'double');
% twoA = zeros(size(F,1),1,'single');
% F_areas = zeros(size(F,1),1,'single');



axis square;graphlims = [-1.1 1.1]; xlim(graphlims);ylim(graphlims); zlim(graphlims);
[R, shear , area_res] = objective(Xo, flag, F,F_areas_relative, DM, invDM, F_areas_o,nvert,[t1 p1], AreaCalcEqual, DMflag);
% [c, ceq] = constraints(Xo, 1)
    options = optimset('MaxFunEvals', 400000,'DiffMaxChange', 1e-2,'DiffMinChange', 1e-8,...
            'DerivativeCheck','off','GradObj','off','GradConstr','off',...
            'MaxIter', maxiter,'Display', 'iter','Diagnostics', 'on','LevenbergMarquardt','on',...
            'TolCon', 1e-6,'TolFun', 1e-16,'TolX', 1e-16, 'MaxSQPIter',30,...
            'LargeScale','on', 'JacobPattern',JacPat,'LineSearchType','quadcubic',...
            'PrecondBandWidth',0,'TypicalX',0 * ones(length(Xo),1),...
            'MaxPCGIter',10);
%             'MaxPCGIter','max(1,floor(numberOfVariables/2))');
% [Xf, resnorm,residual,exitflag,output,lambda, Jacobian] = lsqnonlin(@objective,Xo,[],[],options,1);


[Xf] = lsqnonlin(@objective,Xo,[],[],options,flag, F,F_areas_relative, DM, invDM, F_areas_o,...
    nvert,[t1 p1], AreaCalcEqual, DMflag);
%Xf = fminsearch(@objective2,Xo,options,100, F,F_areas_relative, DM, invDM,
%F_areas_o,nvert);
%[Xf,fval] = fmincon(@objective, Xo, [],[],[],[],[],[],@constraints,options,1);
[R, shear, area_res] = objective(Xf, flag, F,F_areas_relative, DM, invDM, F_areas_o,nvert,[t1 p1], AreaCalcEqual, DMflag);

t = [Xf(1:(nvert))]; p = [Xf(nvert+1:end)];
%t = [Xfixed(1,1); t; Xfixed(2,1)];p = [Xfixed(1,2); p; Xfixed(2,2)];

