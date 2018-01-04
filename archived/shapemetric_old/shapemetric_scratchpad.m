% Write main compoent to select a pixel and update neighbouring pixels
% code is self-contained
% ----------------

% initialize and load test image
close all;
clear all;
im = im2bw(imread('test2.png'));
imshow(im, [], 'InitialMagnification', 'fit');

% extract center and plot
regions = regionprops(im);
centroid = regions(1).Centroid;
center = int8(centroid);
hold on; plot(center(1), center(2), 'rs'); hold off;


% display image data
disp('---- Image_data ----')
disp(im);


% initialise dist_data and display
dist_data = zeros(size(im), 'int8');

disp('---- dist_data ----');
disp(dist_data);

% select pixel from coords and update dist_data
selected_pixel = [center(2), center(1)];
display(selected_pixel);
dist_data = update_dist_old(im, dist_data, selected_pixel);

disp('---- update 1 dist_data ----');
disp(dist_data);

% select new pixel and update again
selected_pixel = selected_pixel + int8([-1,0]);
display(selected_pixel);
dist_data = update_dist_old(im, dist_data, selected_pixel);

dist_data(center(2), center(1)) = 0; %hotfix:relpace center with 0
disp('---- update 2 dist_data ----');
disp(dist_data);



% Next task: write component for selecting pixels



function [dist_data_new] = update_dist_old(im, dist_data_old, selected_pixel)
% Analyse nbhd around selected_pixel in im and add new dist_data

% convert coordinates to index
selected_index = (selected_pixel(2)-1)*size(im,2) + selected_pixel(1);
selected_value = dist_data_old(selected_index);

% initialize output
dist_data_new = dist_data_old;

% loop over rows above and below selected_pixel
x_nbhd = -1:1;
for i = x_nbhd
    % loop over columns left and right of celecter_pixel
    % (excluding selected_pixel)
    if i ~= 0
        y_nbhd = -1:1;
    else
        y_nbhd = -1:2:1;
    end
    
    for j = y_nbhd
        
        % select target pixel
        target_pixel = selected_pixel + int8([i, j]);
        % convert target coordinates to target index
        target_index = (target_pixel(2)-1)*size(im,2)+target_pixel(1);
        
        % if target pixel is outside shape, then set target to -9
        if im(target_index) == 0
            dist_data_new(target_index) = -9;
            
            % otherwise if the target pixel is in the shape and has not been
            % filled in yet, then set target in data_dist to new value
        elseif (im(target_index) ~= 0) && (dist_data_old(target_index) == 0)
            dist_data_new(target_index) = selected_value + 1;
            
        end
    end
end

end