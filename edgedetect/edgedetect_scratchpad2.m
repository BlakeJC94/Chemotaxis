clear all;
close all;

im = imread('shapessm.jpg');
I1 = im;
% figure(1); imshow(I1); title('Original image'); pause(0.5);


I2 = edge(I1, 'canny');
I2 = bwareaopen(I2, 20);
% figure(2); imshow(I2); title('Edge output'); pause(0.5);

I2 = imfill(I2, 'holes');
% figure(2); imshow(I2); title('Imfill output'); pause(0.5);

dilationFactor = 5;
se90 = strel('line', dilationFactor, 90);
se0 = strel('line', dilationFactor, 0);
I3 = imdilate(I2, [se90 se0]);
% figure(3); imshow(I3); title('Imdilate output'); pause(0.5);

I3 = imfill(I3, 'holes');
I3 = bwareaopen(I3, 400);
% figure(3); imshow(I3); title('Imdilate imfill output'); pause(0.5);

seD = strel('diamond',1);
I3 = imerode(I3, seD);
I3 = imerode(I3, seD);
cellOutline = bwperim(I3);
% figure(3); imshow(I3); title('Imerode output'); pause(0.5);


I4 = bwperim(I3);
regions = regionprops(I4);
% centroids = zeros(10,2);
% for i = 1:length(regions)
%     centroids(i,:) = regions(i).Centroid;
% end
% centroids((sum(centroids')>0)',:);

i = 2;
obj_x = floor(regions(i).BoundingBox(1));
obj_y = floor(regions(i).BoundingBox(2));
obj_width = regions(i).BoundingBox(3);
obj_height = regions(i).BoundingBox(4);
obj_centroid = regions(i).Centroid;

obj = I4(obj_y-5:obj_y+obj_height+5, obj_x-5:obj_x+obj_width+5);
imshow(obj);

%https://uk.mathworks.com/matlabcentral/answers/55501-how-get-centroid-contour-distance-at-every-10-degrees
angles = 0:30:360;
x = cosd(angles);
y = sind(angles);

figure;
plot(x,y,'ro');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on;
xCenter = 0;
yCenter = 0;

for k = 1:length(angles)
    
    % Draw line from center to edge of circle
    line([xCenter x(k)], [yCenter y(k)], 'LineWidth', 3);
    
    % calculate angle of line
    angle = atand(y(k)/x(k));
    
    % Convert to degrees
	if y(k) >= yCenter && x(k) >= xCenter
		quadrant = 1;
	elseif y(k) >= yCenter && x(k) <= xCenter
		angle = 180 + angle;
		quadrant = 2;
	elseif y(k) <= yCenter && x(k) < xCenter
		angle = 180 + angle;
		quadrant = 3;
	elseif y(k) <= yCenter && x(k) >= xCenter
		angle = 360 + angle;
		quadrant = 4;
    end
    
    
    
end
