function [S_rot] = translate_rotate_SHP(S, parameters)
%% translate and rotate the group of shapes around their center of mass
%% angle are the Euler angles in the z-x-z convention
%  see:
% http://en.wikipedia.org/wiki/Euler_angles#Euler_angles_as_composition_of_Euler_rotations
% for definition
shift = parameters(1:3);
angles = parameters(4:end);


%% center of mass
[all_F,all_X, cm, r] = get_shps_cms(S);     % obtain center of mass
X = cell2mat(all_X');   % concatenate all X
cm = sum(X,1)./length(X);cm = cm(:)';
% generate the rotation matrix
a = angles(1);
b = angles(2);
g = angles(3);
R = [(cosd(a)*cosd(g)-sind(a)*cosd(b)*sind(g)) ...
               (-cosd(a)*sind(g)-sind(a)*cosd(b)*cosd(g)) ...
                           sind(b)*sind(a);...
     (sind(a)*cosd(g)+cosd(a)*cosd(b)*sind(g)) ...
               (-sind(a)*sind(g)+cosd(a)*cosd(b)*cosd(g)) ...
                           -sind(b)*cosd(a);...
      sind(b)*sind(g)...
               sind(b)*cosd(g)...
                            cosd(b)];
%% translate the center of mass to the origin
[S] = shp_translate2(S, -cm);
% apply the rotation
S_rot = S;
for ix = 1:size(S, 1),
    [xc  yc zc] = get_xyz_clks(S(ix,:));
    C = [xc(:)';yc(:)';zc(:)']';
    cr = [C*R]; tx = [cr(:,1)];ty = [cr(:,2)];tz = [cr(:,3)];
    S_rot(ix,:) = [tx(:)' ty(:)' tz(:)' ];
end
% translate back the center of mass to where is was and add the shift
shift = shift(:)';
[S_rot] = shp_translate2(S_rot, cm - shift);
