% initialize and load test image
close all;
clear all;
im = im2bw(imread('test2.png'));
% imshow(im, [], 'InitialMagnification', 'fit');

% extract center and plot
regions = regionprops(im);
centroid = regions(1).Centroid;
center = round(centroid);
% hold on; plot(center(1), center(2), 'rs'); hold off;


% display image data
% disp('---- Image_data ----')
% disp(im);


% initialise dist_data and display
dist_data = zeros(size(im), 'int8');
%
% disp('---- dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);


% select pixel from coords and update dist_data
selected_pixel = [center(2), center(1)];
selected_index = pix2ind(selected_pixel,im);

display(selected_pixel);
dist_data = update_dist(im, dist_data, selected_index);

old_input = selected_index;
new_input = update_indices(old_input, dist_data);

disp('---- update 1 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);
display(new_input);




for k = 1:length(new_input)
    selected_index = new_input(k);
    selected_pixel = ind2pix(selected_index,im);
    dist_data = update_dist(im, dist_data, selected_index);
    %     disp('---- update 2 dist_data ----');
    %     out = [0:14; (1:14)' dist_data];
    %     disp(out);
end

old_input = new_input;
new_input = update_indices(old_input, dist_data);

disp('---- update 2 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);
display(new_input);




for k = 1:length(new_input)
    selected_index = new_input(k);
    selected_pixel = ind2pix(selected_index,im);
    dist_data = update_dist(im, dist_data, selected_index);
    %     disp('---- update 2 dist_data ----');
    %     out = [0:14; (1:14)' dist_data];
    %     disp(out);
end

old_input = new_input;
new_input = update_indices(old_input, dist_data);

disp('---- update 3 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);
display(new_input);


for k = 1:length(new_input)
    selected_index = new_input(k);
    selected_pixel = ind2pix(selected_index,im);
    dist_data = update_dist(im, dist_data, selected_index);
    %     disp('---- update 2 dist_data ----');
    %     out = [0:14; (1:14)' dist_data];
    %     disp(out);
end

dist_data(center(2), center(1)) = 0; %hotfix:relpace center with 0
old_input = new_input;
new_input = update_indices(old_input, dist_data);

disp('---- update 4 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);
display(new_input);



