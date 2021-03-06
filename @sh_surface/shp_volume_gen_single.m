function [I, block, in] = shp_volume_gen_single(obj, dim, nico, Y_LK, C)
%%%  generate the volume with the object in it
%%%  to generate Y_LK call get_mesh beforehand once with:
%%%  [XF X C Y_LK] = get_mesh(obj, nico)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ==2, nico = 3; [XF X C] = obj.get_mesh(nico);end  % slow, needs to build basis
if nargin ==5, [XF X C] = obj.get_mesh(nico, Y_LK, C);end     % faster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minx = floor(min(X(:,1)));maxx = ceil((max(X(:,1))));
miny = floor(min(X(:,2)));maxy = ceil(max(X(:,2)));
minz = floor(min(X(:,3)));maxz = ceil(max(X(:,3)));
if minx <1,minx = 1;end;if maxx>dim(1), maxx = dim(1);end
if miny <1,miny = 1;end;if maxy>dim(2), maxy = dim(2);end
if minz <1,minz = 1;end;if maxz>dim(3), maxz = dim(3);end

[xi yi zi] = meshgrid(minx:maxx, miny:maxy, minz:maxz);
in = inhull([xi(:) yi(:) zi(:)], X, C);     %%%%------ SLOW -------
in = reshape(in,size(xi));
in = permute(in,[2 1 3]);
block.minx = minx;
block.miny = miny;
block.minz = minz;

block.maxx = maxx;
block.maxy = maxy;
block.maxz = maxz;
if isempty(in),error('empty object');end;

if nargout==3
I = zeros(dim(1),dim(2),dim(3), 'double');
%I(miny:maxy,minx:maxx,minz:maxz) =  I(miny:maxy,minx:maxx,minz:maxz) + double(in);
 I(minx:maxx,miny:maxy,minz:maxz) =  I(minx:maxx,miny:maxy,minz:maxz) + double(in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
