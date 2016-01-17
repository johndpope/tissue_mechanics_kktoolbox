function I = image_pad_z(I, zdimadd)

[xd yd zd] = size(I);
Iz = zeros(xd,yd*zdimadd/2);
I = reshape(I,xd, yd*zd);
I = [Iz I Iz];
I = reshape(I,xd, yd, zd+zdimadd);