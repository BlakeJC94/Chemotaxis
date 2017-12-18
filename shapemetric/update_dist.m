function [dist_data_new] = update_dist(im, dist_data_old, selected_index)
% Analyse nbhd around selected_pixel in im and add new dist_data

% convert index to pixel
selected_pixel = ind2pix(selected_index,im);
selected_value = dist_data_old(selected_index);

% initialize output
dist_data_new = dist_data_old;

if selected_value >= 0
    
    % loop over rows above and below selected_pixel
    i_nbhd = -1:1;
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
            %         display(target_pixel);
            
            % convert target coordinates to target index
            target_index = pix2ind(target_pixel,im);
            
            % if target pixel is outside shape, then set target to -9
            if im(target_index) == 0
                dist_data_new(target_index) = -9;
                
                % otherwise if the target pixel is in the shape and has not been
                % filled in yet, then set target in data_dist to new value
            elseif (im(target_index) ~= 0) && (dist_data_old(target_index) == 0)
                dist_data_new(target_index) = selected_value + 1;
                
            end
            
            
            %         display(dist_data_new);
            
        end
    end
end

end