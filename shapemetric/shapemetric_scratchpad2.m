% initialize and load test image
close all;
clear all;
im = im2bw(imread('test2.png'));
imshow(im, [], 'InitialMagnification', 'fit');

% extract center and plot
regions = regionprops(im);
centroid = regions(1).Centroid;
center = round(centroid);
hold on; plot(center(1), center(2), 'rs'); hold off;


% display image data
disp('---- Image_data ----')
disp(im);


% initialise dist_data and display
dist_data = zeros(size(im), 'int8');

disp('---- dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);


% select pixel from coords and update dist_data
selected_pixel = [center(2), center(1)];
selected_index = pix2ind(selected_pixel,im);

display(selected_pixel);
dist_data = update_dist(im, dist_data, selected_index);

disp('---- update 1 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);



% loop over rows above and below selected_pixel
i_nbhd = -1:1;
k = 0;
A = zeros(1,8);
for i = i_nbhd
    % loop over columns left and right of selected_pixel
    % (excluding selected_pixel)
    if i ~= 0
        j_nbhd = -1:1;
    else
        j_nbhd = -1:2:1;
    end
    for j = j_nbhd
        
        
        % select target pixel
        target_pixel = selected_pixel + [i, j];
        % convert target coordinates to target index
        target_index = pix2ind(target_pixel, im);
        
        k = k+1;
        A(k) = target_index;
        
    end
    
end
    
disp(A) 

for k = 1:8
    selected_index = A(k);
    selected_pixel = ind2pix(selected_index,im)
    dist_data = update_dist(im, dist_data, selected_index);
%     disp('---- update 2 dist_data ----');
%     out = [0:14; (1:14)' dist_data];
%     disp(out);
end

dist_data(center(2), center(1)) = 0; %hotfix:relpace center with 0
disp('---- update 2 dist_data ----');
out = [0:14; (1:14)' dist_data];
disp(out);


