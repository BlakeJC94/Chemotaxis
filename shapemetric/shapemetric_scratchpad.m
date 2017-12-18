close all;
clear all;

im = im2bw(imread('test2.png'));
imshow(im, [], 'InitialMagnification', 'fit');

regions = regionprops(im);
centroid = regions(1).Centroid;

hold on; plot(centroid(1), centroid(2), 'r+'); hold off;

dist_data = zeros(size(im), 'int8');

center = int8(centroid);

selected_pixel = center;
for i = -1:1
    for j = -1:1
        target_index = (selected_pixel(1)+i-1)*size(im,2)+(selected_pixel(2)+j-1);
        
        if im(target_index) == 0
            dist_data(target_index) = -1;
            
        else
            dist_data(target_index) = 1;
            disp('trace');
        end
%         dist_data(selected_pixel(1) + i, selected_pixel(2) + j) = 1;
    end
end