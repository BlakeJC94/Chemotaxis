close all;
clear all;

im = im2bw(imread('test.png'));
imshow(im, [], 'InitialMagnification', 'fit');

regions = regionprops(im);
centroid = regions(1).Centroid;

hold on; plot(centroid(1), centroid(2), 'r+'); hold off;
