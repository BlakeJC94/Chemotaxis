function [dist_data_new input_new] = update_dist_inds(im, dist_data_old, input_old)
% Performs one step of the shapedist transform
% Starting from an initial point 0 inside the shape, look at nbhd (8) and
% update points with +1 if inside the shape and -1 if outside shape.
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
% Blake Cook 2017-12-18



% initialize output
dist_data_new = dist_data_old;
input_new = zeros(1,8*length(input_old));
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
%debug: stop if there are no new indices to update
input_new = sort(input_new(input_new>0));
if isempty(input_new) == 1
    warning('Empty selection');
end

end