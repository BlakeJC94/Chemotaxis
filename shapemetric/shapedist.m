function im_sd = shapedist(im)
% Performs shapedist transform on binary image of single region. Useful for
% finding lengths of shortest paths contained in shape from a center point
% to perimeter points.
%
% Starting from an initial point 0 inside the shape, look at nbhd (8) and
% update points with +1 if inside the shape and -1 if outside shape.
%
% Function will accept binary image of single region and return shapedist
% transfrom integer array
%
% INPUTS:
%   im - binary image of a region
%
% OUTPUTS:
%   im - integer array of shapedist transform
%
% TODO: Make selecting pixels in loop function to analyse more efficient.
%
% Blake Cook 2017-12-19
% ------------------------

%get region data
regions = regionprops(im);
%if there is more than one region, select region with largest area
if length(regions)>1
    disp('More than one region detected, selection largest connected component');
    im_label = labelmatrix(bwconncomp(im));
    
    regions_cell = struct2cell(regions);
    areas = cell2mat(regions_cell(1,:));
    area_index = find(areas == max(areas));
    
    im = ismember(im_label, area_index) > 0;
    regions = regionprops(im);
end

%get center pixel
size_im = size(im);
center_pixel = round(regions(1).Centroid);

%select new center if centroid is not inside the shape
if im(center_pixel(1), center_pixel(2)) == 0
    
    disp('Centroid not located inside shape, finding initial starting point');
    imshow(im,[],'InitialMagnification', 'fit');
    hold on; plot(center_pixel(1), center_pixel(2), 'r+'); hold off;
    
    boundary = bwboundaries(im); %isolate object bondary in image
    boundary = cell2mat(boundary); %convert 1s on image into coordinates
    x = boundary(:,2);
    y = boundary(:,1);
    
    xCenter = center_pixel(1);
    yCenter = center_pixel(2);
    
    [angles, distances] = cart2pol(x-xCenter,y-yCenter);
    
    %get index of boundary point with shortest dist to centroid
    nearest_index = find(distances == min(distances));
    
    %get index of point with smallest angular deviation
    search_angle = angles(nearest_index);
    sorted_angles = sort(abs(angles-search_angle));
    
    %find possible indidies, sort possible distances by length
    possible_index = abs(angles-search_angle) == sorted_angles(2);
    sorted_distances = sort(distances(possible_index));
    
    %get index of point with next smallest distance from centroid
    ray_index = find(distances == sorted_distances(2) & possible_index);
    
    
    %propose new center
    center_pixel(2) = round((x(ray_index) + x(nearest_index))/2);
    center_pixel(1) = round((y(ray_index) + y(nearest_index))/2);
    
%     % -- debug --
%     hold on; plot(center_pixel(2), center_pixel(1), 'm+'); hold off;
%     pause(0.5);
%     
%     % -- debug --
%     plot(x(:), size(im,2) - y(:));
%     hold on;
%     plot(xCenter, size(im,2) - yCenter, 'r+');
%     plot(center_pixel(1), size(im,2) - center_pixel(2), 'm+');
%     hold off;
    
end

%initialise output and input for loop
im_sd = zeros(size(im), 'int8'); %change int size here if image too big)
input = [center_pixel(1), center_pixel(2)];

%set loop maximum to be largest distance of pixels between im boundaries
%plus 2 (include one extra loop for center pixel and another for safety)
% loop_max = max([center_pixel, abs(center_pixel-size(im))+1]) + 1;
loop_max = max(size_im)^2;

for loop = 1:loop_max
    
    [im_sd, new_input, stop_flag] = shapedist_loop(im, im_sd, input);
    input = new_input;
    
%     % -- debug --
%     imshow(mat2gray(im_sd), [], 'InitialMagnification', 'fit');
%     pause(0.5);
    
    if (stop_flag == 1)
        break;
    end
    
end

if loop == loop_max
    warning('reached max loop number, transform may be incomplete');
end

end



function [dist_data_new, input_new, stop_flag] = shapedist_loop(im, dist_data_old, input)
% Performs one step of the shapedist transform.
%
% For a selected point with value n, non-written neighbouring points will
% be updated with vaue n+1. Updated points to n+1 are saved in input_new
%
% Ignores selecting exterior points (value -1)
%
% Function will cycle through selected points in input, look at their
% nbhds and update them as above
%
% INPUTS:
%   im - binary image of a region
%   dis_data_old - int8 array (dims of im) previous step of shapedist
%   input - n x 2 array (double) of pixels to select and update
%
% OUTPUTS:
%   dist_data_new - int8 array next step, after updating input points
%   input_new - m x 2 array (double) of pixels to select on next step
%
% ------------------------


% initialize output
dist_data_new = dist_data_old;
input_new = zeros(8*size(input,1),2);
stop_flag = 0; %signal to stop loop in shapedist
save_index = 0; %row of input_new to update

%for each row in input
for selected_index = 1:size(input,1)
    
    %select pixel and get value
    selected_pixel = input(selected_index,:);
    selected_value = dist_data_old(selected_pixel(1), selected_pixel(2));
    
    %warn if values are too large for int8
    if selected_value >= 127
        warning('Image too big for int8 array, initialise larger array')
    end
    
    %only scan selected pixels with value >= 0 (ignore exterior)
    if selected_value >= 0
        
        %target nbhd pixels
        %loop over rows above and below selected_pixel
        i_nbhd = -1:1;
        for i = i_nbhd
            %loop over columns left and right of selected_pixel
            %(exclude target selected_pixel)
            if i ~= 0
                j_nbhd = -1:1;
            else
                j_nbhd = -1:2:1;
            end
            for j = j_nbhd
                
                %get target_pixel
                target_pixel = selected_pixel + [i, j];
                
                %if target_pixel is outside shape, set target value to -1
                %otherwise, if inside shape and has not been written to yet
                %then update to new value and add to list input_new
                if im(target_pixel(1), target_pixel(2)) == 0
                    dist_data_new(target_pixel(1), target_pixel(2)) = -1;
                    
                elseif (im(target_pixel(1), target_pixel(2)) ~= 0) ...
                        && (dist_data_new(target_pixel(1), target_pixel(2)) == 0)
                    
                    dist_data_new(target_pixel(1), target_pixel(2)) = selected_value + 1;
                    
                    %check if target pixel is on boundary of image
                    boundary_flag = (target_pixel(1) == [1, size(im,1)])...
                        |(target_pixel(2) == [1, size(im,2)]);
                    
                    %if target pixel is not on boundary of image, then add
                    %to list of pixels to update on next loop
                    if boundary_flag ~= 1
                        save_index = save_index+1;
                        input_new(save_index,:) = target_pixel;
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

%remove empty entries from input_new
input_new = input_new(1:save_index,:);
%if empty, stop loop in shapedist
if isempty(input_new) == 1
    stop_flag = 1;
end

end