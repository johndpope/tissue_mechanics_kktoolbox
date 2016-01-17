function [t p bm] = conformal_spherical_mapping_Gu(X,F, t, p)
% % % %% Analyze surface and repair if necessary
bm = [];
disp('Analyzing surface... ');
[X, F] = meshcheckrepair(X,F,'duplicated');
[X, F] = meshcheckrepair(X,F,'isolated');
[X, F] = meshcheckrepair(X,F,'deep', '99');
facecell= finddisconnsurf(F);        % find disconnected surfaces
if size(facecell,2)>2, error('mesh is not suitable for spherical mapping');end
disp('Done!');
%%
tr = TriRep(F,X(:,1), X(:,2), X(:,3));
tr=ConditionMesh(tr);
tri = tr.Triangulation;
F=cell2mat(edgeAttachments(tr,edges(tr))); %% if we don't get an error here then the edges are ok.
%% configure
dt=1E-2; % time step
dE=1E-7; % convergence criterion
Nmax=10000; % maximum number of iterations
% % %% Compute the barycentric map (aka Tutte map)
if nargin>2,
    [x y z] = kk_sph2cart(t,p,1);
    bm=GDSP(tr,[x(:) y(:) z(:)],'bm',{dt dE Nmax}); % [] means no initialization was provided
else
     bm=GDSP(tr,[],'bm',{dt dE Nmax}); % [] means no initialization was provided
end
%% Visualize the barycentric map
figure('color','w')
ha=axes;
h=trimesh(TriRep(tri,bm));
set(h,'EdgeColor','b'), axis off, axis equal, drawnow
set(get(ha,'Title'),'String','Barycentric Map',...
'FontSize',35,'FontWeight','bold')

%% Compute the conformal map (for genus-0 surfaces harmonic and conformal maps are equivalent ).
if nargin>2,
    [x y z] = kk_sph2cart(t,p,1);
    cm=GDSP(tr,[x(:) y(:) z(:)],'cm',{dt dE Nmax}); % use input to initialize the search
else
    cm=GDSP(tr,bm,'cm',{dt dE Nmax}); %
end
%% Visualize the conformal map
figure('color','w')
ha=axes;
h=trimesh(TriRep(tri,cm));
set(h,'EdgeColor','b'), axis off, axis equal, drawnow
set(get(ha,'Title'),'String','Conformal Map',...
'FontSize',35,'FontWeight','bold');

%% Visualize the distribution of angle distortions
%AngleDistortionPlot(tr,bm,'Barycentic Map Ang Dist Err');
AngleDistortionPlot(tr,cm,'Conformal Map Ang Dist Err');
%% output
[t p r] = kk_cart2sph(cm(:,1), cm(:,2), cm(:,3));






