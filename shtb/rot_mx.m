function R = rot_mx(ang,flag)
%% return a cosine rotation matrix angles in radians
ca = cos(ang);
sa = sin(ang);
if flag == 1,    %% rotate around z
    R = [ca -sa 0;sa ca 0;0 0 1];
elseif flag == 2,       %% rotate around y
    R = [ca 0 sa;0 1 0;-sa 0 ca];
elseif flag == 3,       %% rotate around x
    R = [1 0 0;0 ca -sa;0 sa ca];
end