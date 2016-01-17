function h = overlay_image_segmentation2(I, all_X, all_F)
printing = 1;
%%% overlay image I and segmentation
h = figure;
im = mat2gray(sum(I,3));
figure(1);imshow(squeeze(im));axis equal
xlim([0 size(im,2)]);ylim([0 size(im,1)]);
if printing, print -dtiff -r600 projection_view_1.tif;end
figure(2);imshow(im);axis equal
xlim([0 size(im,2)]);ylim([0 size(im,1)]);

%%
all_X1 = all_X;
for ix = 1:length(all_X)
    xx = all_X{ix};
    xx(:,3) = 0;
    all_X1{ix} = xx;
end
%%
hold on;plot_all_X(all_X1,all_F);
if printing, print -dtiff -r600 overlay_view_1.tif;end
%% look from different angle
im = mat2gray(sum(I,1));
figure(3);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,2)]);
if printing, print -dtiff -r600 projection_view_2.tif;end

figure(4);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,2)]);

%%
all_X2 = all_X;
for ix = 1:length(all_X)
    xx = all_X{ix};
    temp = xx(:,1);
    xx(:,1) = xx(:,3);  % x becomes the new y
    xx(:,2) = temp;  % switch z and y
    xx(:,3) = 0;
    all_X2{ix} = xx;
end
%%
hold on;plot_all_X(all_X2,all_F);
if printing, print -dtiff -r600 overlay_view_2.tif;end
%% yet from another different angle
im = mat2gray(sum(I,2));
figure(5);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,1)]);
if printing, print -dtiff -r600 projection_view_3.tif;end
figure(6);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,1)]);

%%
all_X3 = all_X;
for ix = 1:length(all_X)
    xx = all_X{ix};
    temp = xx(:,2);
    xx(:,1) = xx(:,3);  % x becomes the new y
    xx(:,2) = temp;  % switch z and y
    xx(:,3) = 0;
    all_X3{ix} = xx;
end
%%
hold on;plot_all_X(all_X3,all_F);
if printing, print -dtiff -r600 overlay_view_3.tif;end