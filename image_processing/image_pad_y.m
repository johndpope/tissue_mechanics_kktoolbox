function I = image_pad_y(I, ydimadd)

%%% change the axis, pad along z and then change them back
I = permute(I,[1 3 2]);
I = image_pad_z(I,ydimadd);
I = permute(I,[1 3 2]);