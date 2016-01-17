function [all_X_new, cm] = translate_rotate_XF(all_X, parameters)
%% translate and rotate the group of shapes around their center of mass
%% angle are the Euler angles in the z-x-z convention
%  see:
% http://en.wikipedia.org/wiki/Euler_angles#Euler_angles_as_composition_of_Euler_rotations
% for definition
shift = parameters(1:3);
angles = parameters(4:end);
convert_to_cell = 0;
if ~iscell(all_X), convert_to_cell =1;X{1} = all_X; all_X = X;end
X = cell2mat(all_X');   % concatenate all X
%% center of mass
cm = sum(X,1)./length(X);cm = cm(:)';
% generate the rotation matrix
a = angles(1);
b = angles(3);
g = angles(1);
R = [(cosd(a)*cosd(g)-sind(a)*cosd(b)*sind(g)) ...
               (-cosd(a)*sind(g)-sind(a)*cosd(b)*cosd(g)) ...
                           sind(b)*sind(a);...
     (sind(a)*cosd(g)+cosd(a)*cosd(b)*sind(g)) ...
               (-sind(a)*sind(g)+cosd(a)*cosd(b)*cosd(g)) ...
                           -sind(b)*cosd(a);...
      sind(b)*sind(g)...
               sind(b)*cosd(g)...
                            cosd(b)];
% translate the center of mass to the origin
X = X-cm(ones(1,size(X,1)), :);
% apply the rotation
Xp = (X*R);
% translate back the center of mass to where is was and add the shift
shift = shift(:)';
Xp = Xp+cm(ones(1,size(X,1)), :)+ shift(ones(1,size(X,1)), :);
% generate the new set of objects
counter = 1;
all_X_new = all_X;
for ix = 1:size(all_X,2),
    l = size(all_X{ix},1);
    all_X_new{ix} = Xp(counter:l+counter-1, 1:3);
    counter = counter + l;
%     Xp(1:l, :) = [];
end

if convert_to_cell, all_X_new = cell2mat(all_X_new);end