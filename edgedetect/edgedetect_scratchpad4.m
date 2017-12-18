clear all;
close all;

%% Import test imag, filter and output label matrix for connected regions
im = imread('shapessm.jpg');
I1 = im;
figure(1); imshow(I1); title('Original image'); pause(0.5);

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
% figure(3); imshow(I3); title('Imerode output'); pause(0.5);

cc = bwconncomp(I3);
labeled = labelmatrix(cc);
figure(1); imshow(label2rgb(labeled)); title('rgblabelled regions'); pause(0.5);


%% Extract region properties
regions = regionprops(I3);
i = 2; % Looking at floppy disk
obj = ismember(labeled, i) > 0;
figure(1); imshow(obj); title(['Object i = ' num2str(i)]); pause(0.5);

obj_boundary = bwboundaries(obj); %isolate object bondary in image
obj_boundary = cell2mat(obj_boundary); %convert 1s on image into coordinates
x = obj_boundary(:,2);
y = obj_boundary(:,1);

xCenter = regions(i).Centroid(1);
yCenter = regions(i).Centroid(2);

[angles, distances] = cart2pol(x-xCenter,y-yCenter);
angles = angles*180/pi + 180;
[angles, index] = sort(angles);
distances = distances(index);

% PLOT DIST EVERY angle_step DEGREES
angle_step = 1;
if mod(360,angle_step) ~= 0
    warning('angle_step not an integer divisor of 360');
end

angles_vec = 0:angle_step:359.9999;
dist_vec = zeros(size(angles_vec));

angles = floor(angles / angle_step) * angle_step;
for k = 1:length(angles_vec)
    dist_vec(k) = distances(find(angles >= angles_vec(k), 1));
end
figure(2); subplot(1,2,1);
plot(angles_vec,dist_vec, 'bo-');
xlabel('angles (degrees)');
ylabel('dist');

subplot(1,2,2);
polarplot(angles_vec*(pi/180),dist_vec, 'bo-');

% -------------

% distances = sqrt((x - xCenter).^2 + (y - yCenter).^2);
% angles = zeros(size(x));
% for k = 1 : length(x)
% 	angles(k) = atand((y(k)-yCenter) / (x(k)-xCenter));
% 
% 	if y(k) >= yCenter && x(k) >= xCenter
% 		quadrant = 1;
% 	elseif y(k) >= yCenter && x(k) <= xCenter
% 		angles(k) = 180 + angles(k);
% 		quadrant = 2;
% 	elseif y(k) <= yCenter && x(k) < xCenter
% 		angles(k) = 180 + angles(k);
% 		quadrant = 3;
% 	elseif y(k) <= yCenter && x(k) >= xCenter
% 		angles(k) = 360 + angles(k);
% 		quadrant = 4;
% 	end
% end
% [angles, index] = sort(angles);
% distances = distances(index);
