function [allX allF] = write_shps_2_amira(S,gdim, fn)
%%% generate a wave front .obj file from a group of spherical harmonics
%%% shapes 
allX = [];
allF = [];
verbose = 0;
if size(S,2) ==1,S = S(:)';end;% check whether S has a dimension which is one. i.e. we have only one shape
L_max = round(sqrt(size(S,2)/3)-1);
P = partsphere(gdim^2);x = P(1,:);y = P(2,:);z = P(3,:);
[t p r] = kk_cart2sph(x,y,z);[x y z] = kk_sph2cart(t,p,1);
X = [x(:) y(:) z(:)]; [C] = convhulln(X, {'Qt'});
Y_LK = get_basis(t',p',gdim,max([L_max L_max L_max]));
counter = 0;

mkdir('./temp');
cd temp
for ix = 1:size(S,1), %loop over the shapes
    [xclks yclks zclks] = get_xyz_clks(S(ix,:));
    X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
    allX = [X];
    allF = [C];
    
    % write stl file
    str1 = sprintf('%s_%d.raw',fn, ix);str2 = sprintf('%s_%d.stl',fn, ix);mesh2raw(str1,allX,allF);Raw2Stl(str1, str2);
    delete(str1);
    disp(['Exported stl format: ' num2str(ix) ' of ' num2str(size(S,1))]);
    counter = counter+1;
end

%% write hx file for amira
path = sprintf('./%s.hx', fn);
[fid,Msg] = fopen(path,'wt');
if fid == -1, error(Msg); end
% write header
fprintf(fid,'# Amira Script\nremove -All\n\n# Create viewers\nviewer setVertical 0\n\nviewer 0 setBackgroundMode 1\nviewer 0 setBackgroundColor 0.06 0.13 0.24\nviewer 0 setBackgroundColor2 0.72 0.72 0.78\nviewer 0 setTransparencyType 5\nviewer 0 setAutoRedraw 0\nviewer 0 show\nmainWindow show');

% define the objects to be read
for ix = 1:size(S,1), %loop over the shapes
    str = sprintf('%s_%d.stl',fn, ix);
    fprintf(fid,'\n\nset hideNewModules 0\n[ load -stl ${SCRIPTDIR}/%s ] setLabel %s\n%s setIconPosition 20 10\n%s fire\n%s LevelOfDetail setMinMax -1 -1\n%s LevelOfDetail setButtons 1\n%s LevelOfDetail setIncrement 1\n%s LevelOfDetail setValue -1\n%s LevelOfDetail setSubMinMax -1 -1\n%s fire\n%s setViewerMask 16383\n',str,str, str,str,str, str, str, str);
end

% create the data directory for each object
for ix = 1:size(S,1), %loop over the shapes
    str = sprintf('%s_%d.stl',fn, ix);
    fprintf(fid,'\nset hideNewModules 1\ncreate HxDataDirectory {DataDirectory_%s}\nDataDirectory_%s setIconPosition 0 0\nDataDirectory_%s setInvisible 1\nDataDirectory_%s parentDir connect Data\nDataDirectory_%s attachedData connect %s\nDataDirectory_%s fire\nDataDirectory_%s setViewerMask 16383',...
        str, str, str, str, str, str, str, str);
end


% define the views for each object
for ix = 1:size(S,1), %loop over the shapes
    str = sprintf('%s_%d.stl',fn, ix);
    sv = sprintf('SurfaceView%s',num2str(ix));
    fprintf(fid,'\nset hideNewModules 0\ncreate HxDisplaySurface {%s}\n%s setIconPosition 367 10\n%s data connect %s\n%s colormap setDefaultColor 1 0.1 0.1\n%s colormap setDefaultAlpha 0.500000\n%s fire\n%s drawStyle setValue 1\n%s fire\n%s drawStyle setSpecularLighting 1\n%s drawStyle setTexture 1\n%s drawStyle setAlphaMode 1\n%s drawStyle setNormalBinding 0\n%s drawStyle setSortingMode 1\n%s drawStyle setLineWidth 0.000000\n%s drawStyle setOutlineColor 0 0 0.2\n%s textureWrap setIndex 0 1\n%s cullingMode setValue 0\n%s selectionMode setIndex 0 0\n%s Patch setMinMax 0 1\n%s Patch setButtons 1\n%s Patch setIncrement 1\n%s Patch setValue 0\n%s Patch setSubMinMax 0 1\n%s BoundaryId setIndex 0 -1\n%s materials setIndex 0 1\n%s materials setIndex 1 0\n%s colorMode setIndex 0 0\n%s baseTrans setMinMax 0 1\n%s baseTrans setButtons 0\n%s baseTrans setIncrement 0.1\n%s baseTrans setValue 0.8\n%s baseTrans setSubMinMax 0 1\n%s VRMode setIndex 0 0\n%s fire\n%s hideBox 1\n{%s} selectTriangles zab HIJMPLPPBPDLGAGAGAOAAHAAIKKOBHPI\n%s fire\n%s setViewerMask 16383\n%s setShadowStyle 1\n%s setPickable 1\n',...
        sv,sv, sv,str,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,...
        sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv,sv);
end


% the bottom part of the script is set here.

fprintf(fid,'\nset hideNewModules 0\n\nviewer 0 setCameraOrientation 1 0 0 3.14159\nviewer 0 setCameraPosition 52.9061 49.4087 -69.1548\nviewer 0 setCameraFocalDistance 132.747\nviewer 0 setCameraNearDistance 116.993\nviewer 0 setCameraFarDistance 148.532\nviewer 0 setCameraType perspective\nviewer 0 setCameraHeightAngle 44.9023\nviewer 0 setAutoRedraw 1\nviewer 0 redraw\n');
fclose(fid);

cd ..

