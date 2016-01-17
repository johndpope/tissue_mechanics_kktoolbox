function h = overlay_image_segmentation4(I, all_X, all_F)
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
