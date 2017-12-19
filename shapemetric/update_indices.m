function new_input = update_indices(old_input, dist_data)
% return new selection by scanning nbhds of old selection for update_dist
% input and output are row vectors of indices
% (to be used after update_dist)

new_input = zeros(1,8*length(old_input));

i_nbhd = -1:1;
k = 0;
% loop over old indices and select pixel
for old_index = old_input
    
    selected_pixel = ind2pix(old_index, dist_data);
    
    boundary_flag = (selected_pixel(1) == [1, size(dist_data,1)]) |...
        (selected_pixel(2) == [1, size(dist_data,2)]);
    
    % if selected pixel is not on boundary, continue. Otherwise loop over
    if boundary_flag ~= 1
        
        % loop over rows above and below selected pixel
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
                target_index = pix2ind(target_pixel, dist_data);
                
                % If the targeted pixel has been filled in by a positive
                % integer last time update_dist was used, then add to list of
                % new inputs for update_dist
                if dist_data(target_index) > 0
                    k = k+1;
                    new_input(k) = target_index;
                end
                
            end
        end
    end
end

new_input = sort(new_input(new_input>0));
if isempty(new_input) == 1
    error('Empty selection');
end



end
