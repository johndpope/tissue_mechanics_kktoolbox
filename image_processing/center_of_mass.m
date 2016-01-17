function [av_i, av_j, av_k] = center_of_mass(im)
% calculate the center of mass of a grey value image
% please vectorize
im = squeeze(im);
av_i = 0;
av_j = 0;
av_k = 0;
counter = 0;
total = size(im,1);
h = waitbar(0,'Determining center of mass .....');
for ix = 1:size(im,1),
    counter = counter + 1;
    for jx = 1:size(im,2),
        for kx = 1:size(im,3);
            
            av_i = av_i + ix.*im(ix,jx, kx);
            av_j = av_j + jx.*im(ix,jx, kx);
            av_k = av_k + kx.*im(ix,jx, kx);
            
        end
    end
    waitbar(counter/total)
end
close(h);drawnow;

av_i = av_i/sum(im(:));
av_j = av_j/sum(im(:));
av_k = av_k/sum(im(:));
