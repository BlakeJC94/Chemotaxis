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
% figure(3); imshow(I3); title('Imerode output'); pause(0.5);

cc = bwconncomp(I3);
labeled = labelmatrix(cc);
imshow(label2rgb(labeled));

regions = regionprops(I3);
% centroids = zeros(10,2);
% for i = 1:length(regions)
%     centroids(i,:) = regions(i).Centroid;
% end
% centroids((sum(centroids')>0)',:);

i = 2;
obj = ismember(labeled, i) > 0;

obj_x = floor(regions(i).BoundingBox(1));
obj_y = floor(regions(i).BoundingBox(2));
obj_width = regions(i).BoundingBox(3);
obj_height = regions(i).BoundingBox(4);
obj_centroid = regions(i).Centroid;

%obj = I3(obj_y-5:obj_y+obj_height+5, obj_x-5:obj_x+obj_width+5);
imshow(obj); 

Boundary = bwboundaries(obj);
Boundary = cell2mat(Boundary);

x = Boundary(:,2);
y = Boundary(:,1);

xCenter = obj_centroid(1);
yCenter = obj_centroid(2);


distances = sqrt((x - xCenter).^2 + (y - yCenter).^2);


angles = zeros(size(x));
for k = 1 : length(x)

	% Calculate the angle of that line
	angles(k) = atand((y(k)-yCenter) / (x(k)-xCenter));
	% Convert to 0-360
	if y(k) >= yCenter && x(k) >= xCenter
		quadrant = 1;
	elseif y(k) >= yCenter && x(k) <= xCenter
		angles(k) = 180 + angles(k);
		quadrant = 2;
	elseif y(k) <= yCenter && x(k) < xCenter
		angles(k) = 180 + angles(k);
		quadrant = 3;
	elseif y(k) <= yCenter && x(k) >= xCenter
		angles(k) = 360 + angles(k);
		quadrant = 4;
	end
end

plot(angles,distances, 'b-');
xlabel('angles');
ylabel('dist');

[angles, index] = sort(angles);
distances = distances(index);
plot(angles,distances, 'b-');
xlabel('angles');
ylabel('dist');


%https://uk.mathworks.com/matlabcentral/answers/55501-how-get-centroid-contour-distance-at-every-10-degrees
% angles = 0:30:360;
% x = cosd(angles);
% y = sind(angles);
% 
% figure;
% plot(x,y,'ro');
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% grid on;
% xCenter = 0;
% yCenter = 0;
% 
% for k = 1:length(angles)
%     
%     % Draw line from center to edge of circle
%     line([xCenter x(k)], [yCenter y(k)], 'LineWidth', 3);
%     
%     % calculate angle of line
%     angle = atand(y(k)/x(k));
%     
%     % Convert to degrees
% 	if y(k) >= yCenter && x(k) >= xCenter
% 		quadrant = 1;
% 	elseif y(k) >= yCenter && x(k) <= xCenter
% 		angle = 180 + angle;
% 		quadrant = 2;
% 	elseif y(k) <= yCenter && x(k) < xCenter
% 		angle = 180 + angle;
% 		quadrant = 3;
% 	elseif y(k) <= yCenter && x(k) >= xCenter
% 		angle = 360 + angle;
% 		quadrant = 4;
%     end
%     
% end

% -------

% For each blob, get its boundaries and find the distance from the centroid to each boundary point.
% % 	for k = 1 : numRegions
% % 		% Extract just this blob alone.
% % 		thisBlob = ismember(labeledImage, k) > 0;
% % 		if promptUser % Let user see the image.
% % 			cla;
% % 			imshow(thisBlob);
% % 		end
% % 		% Find the boundaries
% % 		thisBoundary = bwboundaries(thisBlob);
% % 		thisBoundary = cell2mat(thisBoundary); % Convert from cell to double.
% % 		% Get x and y
% % 		x = thisBoundary(:, 2);
% % 		y = thisBoundary(:, 1);
% % 		% Get the centroid
% % 		xCenter = blobMeasurements(k).Centroid(1);
% % 		yCenter = blobMeasurements(k).Centroid(2);
% % 		% Compute distances
% % 		distances = sqrt((x - xCenter).^2 + (y - yCenter).^2);
% % 		if promptUser % Let user see the curves.
% % 			% Plot the distances.
% % 			plot(distances, 'b-', 'LineWidth', 3);
% % 			grid on;
% % 			message = sprintf('Centroid to perimeter distances for shape #%d', k);
% % 			title(message, 'FontSize', 15);
% % 			% Scale y axis
% % 			yl = ylim();
% % 			ylim([0, yl(2)]); % Set lower limit to 0.
% % 		end
% % 		
% % 		% Find the range of the peaks
% % 		peakRange = max(distances) - min(distances);
% % 		minPeakHeight = 0.5 * peakRange;
% % 		% Find the peaks
% % 		[peakValues, peakIndexes] = findpeaks(distances, 'MinPeakProminence', minPeakHeight);
% % 		% Find the valueys.
% % 		[valleyValues, valleyIndexes] = findpeaks(-distances, 'MinPeakProminence', minPeakHeight);
% % 		numVertices(k) = max([length(peakValues), length(valleyValues)]);
% % 		% Circles seem to have a ton of peaks due to the very small range and quanitization of the image.
% % 		% If the number of peaks is more than 10, make it zero to indicate a circle.
% % 		if numVertices(k) > 10
% % 			numVertices(k) = 0;
% % 		end
% % 		
% % 		if promptUser % Let user see the curves.
% % 			% Plot the peaks.
% % 			hold on;
% % 			plot(peakIndexes, distances(peakIndexes), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
% % 			
% % 			% Plot the valleys.
% % 			hold on;
% % 			plot(valleyIndexes, distances(valleyIndexes), 'rv', 'MarkerSize', 10, 'LineWidth', 2);
% % 			
% % 			message = sprintf('Centroid to perimeter distances for shape #%d.  Found %d peaks.', k, numVertices(k));
% % 			title(message, 'FontSize', 20);
% % 			
% % 			% The figure un-maximizes each time when we call cla, so let's maximize it again.
% % 			% Set up figure properties:
% % 			% Enlarge figure to full screen.
% % 			set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % 			% Get rid of tool bar and pulldown menus that are along top of figure.
% % 			set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% % 			% Give a name to the title bar.
% % 			set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
% % 			
% % 			% Let user see this shape's distances plotted before continuing.
% % 			promptMessage = sprintf('Do you want to Continue processing,\nor Cancel processing?');
% % 			titleBarCaption = 'Continue?';
% % 			button = questdlg(promptMessage, titleBarCaption, 'Continue', 'Cancel', 'Continue');
% % 			if strcmpi(button, 'Cancel')
% % 				promptUser = false;
% % 			end
% % 		end
% % 	end

