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
% TODO: Propse new center if centroid is located outside region
% TODO: handle larger images 
% 
% Blake Cook 2017-12-19
% ------------------------

%get region data, stop if there is more than one region
regions = regionprops(im);
if length(regions)>1
    error('More than one region detected');
end

%get center pixel
center_pixel = round(regions(1).Centroid);
center_index = pix2ind(center_pixel,im);

%temp: stop if centroid is not inside shape
if im(center_index) == 0
    error('Centroid not located in shape');
end

%temp: stop if shorted distacnce between center and image boundary larger
%than 127, since im_sd is int8 array
if min([center_pixel-1, abs(center_pixel-size(im))]) > 127
    error('Image dimensions too large for int8 array');
end

%initialise output and input for loop
im_sd = zeros(size(im), 'int32'); %change int size here if image too big)
input = center_index;

%set loop maximum to be largest distance of pixels between im boundaries
%plus 2 (include one extra loop for center pixel and another for safety)
loop_max = max([center_pixel, abs(center_pixel-size(im))+1]) + 1;

for loop = 1:loop_max
    
    [im_sd, new_input, stop_flag] = shapedist_loop(im, im_sd, input);
    input = new_input;
    
%     display(loop)
%     imshow(mat2gray(im_sd), [], 'InitialMagnification', 'fit');

    if (stop_flag == 1)
        break;
    end
    
end

end



function [dist_data_new, input_new, stop_flag] = shapedist_loop(im, dist_data_old, input_old)
% Performs one step of the shapedist transform.
%
% For a selected point with value n, non-written neighbouring points will
% be updated with vaue n+1. Updated points to n+1 are saved in input_new
% 
% Ignores selecting exterior points (value -1)
%
% Function will cycle through selected points in input_old, look at their
% nbhds and update them as above
%
% INPUTS:
%   im - binary image of a region
%   dis_data_old - int8 array (dims of im) previous step of shapedist
%   input_old - 1 x n vector (double) of indices to select and update
%
% OUTPUTS: 
%   dist_data_new - int8 array next step, after updating input_old points
%   input_new - row vector (double) of indices to select on next step
%
% ------------------------


% initialize output
dist_data_new = dist_data_old;
input_new = zeros(1,8*length(input_old));
stop_flag = 0;
tmp = 0;

%for each index input
for selected_index = input_old
    
    %convert index to pixel 
    selected_pixel = ind2pix(selected_index,im);
    selected_value = dist_data_old(selected_index);
    
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
                
                %convert target_pixel to target_index
                target_index = pix2ind(target_pixel,im);
                
                %if target_pixel is outside shape, set target value to -1
                %otherwise, if inside shape and has not been written to yet
                %then update to new value and add to list input_new
                if im(target_index) == 0
                    dist_data_new(target_index) = -1;
                elseif (im(target_index) ~= 0) && (dist_data_new(target_index) == 0)
                    dist_data_new(target_index) = selected_value + 1;
                    
                    %check if target pixel is on boundary of image
                    boundary_flag = (target_pixel(1) == [1, size(im,1)])...
                        |(target_pixel(2) == [1, size(im,2)]);
                    
                    %if target pixel is not on boundary of image, then add
                    %to list of pixels to update on next loop
                    if boundary_flag ~= 1
                        tmp = tmp+1;
                        input_new(tmp) = target_index;
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

%remove empty entries from input_new
%debug: warn if there are no new indices to update
input_new = sort(input_new(input_new>0));
if isempty(input_new) == 1
    stop_flag = 1;
end

end 


function ind = pix2ind(pix, im)
% converts pix (1x2 row,column) into index
ind = (pix(2)-1)*size(im,1) + pix(1);
end

function pix = ind2pix(ind, im)
% converts pix (1x2 row,column) into index
pix = [mod(ind, size(im,1)), ...
    (ind - mod(ind, size(im,1)))/size(im,1) + 1];
end