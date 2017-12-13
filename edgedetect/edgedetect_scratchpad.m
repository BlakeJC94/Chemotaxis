clear all;
close all;

im = imread('shapessm.jpg');
figure(1); imshow(im); title('Original image');

edge_detect = 1;
loop = 0;
cannyfact = 0.3;

while edge_detect == 1
    edge_detect = 0;
    loop = loop + 1;
    
    edges = edge(im, 'canny', cannyfact);
    
%     cc = bwconncomp(edges);
%     labeled = labelmatrix(cc);
%     cc = label2rgb(labeled);
%     figure; imshow(cc); title('Labelled connected components');
    
    edgesfill = imfill(edges, 'holes');
    
%     cc = bwconncomp(edgesfill);
%     labeled = labelmatrix(cc);
%     cc = label2rgb(labeled);
%     figure; imshow(cc); title('Labelled connected components FILL');
    
    regions = regionprops(edgesfill);
    regions_outline = regionprops(bwperim(edgesfill));
    
    if length(regions) ~= length(regions_outline)
        error('number of regions not conserved during filling operation');
    end
    
    for i = 1:length(regions)
        area_fill = regions(i).Area;
        area_outline = regions_outline(i).Area;
        
        if area_outline/area_fill > 0.9
            disp('trace');
            cannyfact = cannyfact*0.9;
            edge_detect = 1;
            
            break;
            %        obj_x = floor(regions_outline(i).BoundingBox(1));
            %        obj_y = floor(regions_outline(i).BoundingBox(2));
            %        obj_width = regions_outline(i).BoundingBox(3);
            %        obj_height = regions_outline(i).BoundingBox(4);
            %
            %        obj = edgesfill(obj_y:obj_y+obj_height, obj_x:obj_x+obj_width);
            %        imshow(obj);
            
        end
        
    end
    
    disp(edge_detect);
    figure(2); imshow(edgesfill); title(['edge detection, canny ' num2str(cannyfact) 'FILL']);
    pause;
    
    if loop > 20
        disp('too many loops')
        break;
    end
end


disp(loop)
