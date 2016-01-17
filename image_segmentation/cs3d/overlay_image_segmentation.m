function h = overlay_image_segmentation(I, all_X, all_F)
%%% overlay image I and segmentation
h = figure;
im = mat2gray(sum(I,3));
subplot(3,2,1);imshow(squeeze(im));axis equal
xlim([0 size(im,2)]);ylim([0 size(im,1)]);
subplot(3,2,2);imshow(im);axis equal
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

%% look from different angle
im = mat2gray(sum(I,1));
subplot(3,2,3);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,2)]);
subplot(3,2,4);imshow(squeeze(im));axis equal
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


%% yet from another different angle
im = mat2gray(sum(I,2));
subplot(3,2,5);imshow(squeeze(im));axis equal
xlim([0 size(im,3)]);ylim([0 size(im,1)]);
subplot(3,2,6);imshow(squeeze(im));axis equal
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