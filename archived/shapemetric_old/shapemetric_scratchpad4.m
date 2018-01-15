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


% initialise dist_data
dist_data = zeros(size(im), 'int8');
%
% disp('---- dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);


% select pixel from coords and update dist_data
selected_pixel = [center(2), center(1)];
display(selected_pixel);

selected_index = pix2ind(selected_pixel,im);

old_input = selected_index;
new_input = 0;


loop = 0;
while (~isempty(new_input)) || (loop>10)
    
    loop = loop + 1;
    [dist_data, new_input] = update_dist_inds(im, dist_data, old_input);
    old_input = new_input;
    
    
    disp(['---- update ' num2str(loop) ' dist_data ----']);
    out = [0:14; (1:14)' dist_data];
    disp(out);
    
end


% [dist_data new_input] = update_dist_inds(im, dist_data, old_input);
% % 
% % old_input = selected_index;
% % new_input = update_indices(old_input, dist_data);
% 
% disp('---- update 1 dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);
% display(new_input);
% 
% 
% 
% old_input = new_input;
% [dist_data, new_input] = update_dist_inds(im, dist_data, old_input);
% 
% disp('---- update 2 dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);
% display(new_input);
% 
% 
% old_input = new_input;
% [dist_data, new_input] = update_dist_inds(im, dist_data, old_input);
% 
% disp('---- update 3 dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);
% display(new_input);
% 
% 
% old_input = new_input;
% [dist_data, new_input] = update_dist_inds(im, dist_data, old_input);
% dist_data(center(2), center(1)) = 0; %hotfix:relpace center with 0
% 
% disp('---- update 4 dist_data ----');
% out = [0:14; (1:14)' dist_data];
% disp(out);
% display(new_input);


