function AnalyseDictyShape
% Function for analysing the images from Dicty data
% Based on code from David M. Richards - 14/08/2017
% Blake Cook - 6/12/2017
% ======================== %

%% Analysis paramters:
%  - CHANNEL  - use channel number CHANNEL(1) out of total of CHANNEL(2) channels
%  - FRAME_RANGE - which frames to analysis; [-1 Inf] => all frames
%  - ROI - region of interest as [x1 y1 x2 y2]; [-1 -1 Inf Inf] => full image
%  - REGIONS_TO_IGNORE - regions to remove from final image as [x1 y1 x2 y2]

addpath('functions/')
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
clf;
set(gcf, 'Position', get(0,'Screensize'));

% *** File "DictyElectrotaxis_171116_001" - 953Mb - 271 frames - 3 channels - 1280x960
DIR = 'data/'; %B: Directory changed to work on my laptop
FILENAME = 'DictyElectrotaxis_171116_001.tif';
CHANNEL = [3 3];
FRAME_RANGE = [155 271]; %B: looking at few frames for now, default "[155 271]"
FRAME_JUMP = 1;
ROI = [340 500 880 800];
REGIONS_TO_IGNORE = [95 5 185 100; 190 1 280 25; 320 1 350 20; 310 140 390 210; 480 225 510 260];
AREA_LIMITS = [500 10000];


%% Plot parameters:
%  - PLOT_NUMS - plot numbers to show
%  - PLOT_MARGIN - margins around each subplot
%  - PLOT_NUM_COLS_LIST, PLOT_NUM_ROWS_LIST - number of rows/rows for subplot
%  - IMAGE_NAMES - names for images

PLOT_NUMS = 1:4;
PLOT_MARGIN_X = 0.003;
PLOT_MARGIN_Y = 0;
PLOT_NUM_ROWS_LIST = [ 1 1 2 2 2 2 3 3 ];
PLOT_NUM_COLS_LIST = [ 1 2 2 2 3 3 3 3 ];
IMAGE_NAMES = {'I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7'};


%% Read header info

file_info = imfinfo([DIR FILENAME]);
im_width = file_info(1).Width;
im_height = file_info(1).Height;
num_images = numel(file_info);

if mod(num_images,CHANNEL(2))~=0
    error(['Total number of frames (' num2str(num_images) ') not multiple of number of channels (' num2str(CHANNEL(2)) ')']);
end
num_frames = num_images/CHANNEL(2);

if (FRAME_RANGE(1)==-1); FRAME_RANGE(1)=1; end
if (FRAME_RANGE(2)==Inf); FRAME_RANGE(2)=num_frames; end

num_frames_range = numel(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2));

if (ROI(1)==-1); ROI(1)=1; end
if (ROI(2)==-1); ROI(2)=1; end
if (ROI(3)==Inf); ROI(3)=im_width; end
if (ROI(4)==Inf); ROI(4)=im_height; end

ROI_width  = ROI(3) - ROI(1) + 1;
ROI_height = ROI(4) - ROI(2) + 1;

num_plots = numel(PLOT_NUMS);


%% Declare variables:
%  - im - holds all images
%  - im_full - holds full frame

im = zeros(ROI_height,ROI_width,num_frames,'uint8');
im_full = zeros(im_height,im_width);


%% Read in images

tif_obj = Tiff([DIR FILENAME],'r');
for i = 1:num_frames
    tif_obj.setDirectory(CHANNEL(1)+CHANNEL(2)*(i-1));
    im_full = tif_obj.read();
    im(:,:,i) = im_full(ROI(2):ROI(4),ROI(1):ROI(3));
end
tif_obj.close();

%% Find background image

% im_bg = zeros(ROI_height, ROI_width, 'double');
% for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
%     for i = 1:ROI_height
%         for j = 1:ROI_width
%             im_bg(i,j) = im_bg(i,j) + double(im(i,j,frameNum));
%         end
%     end
% end
% im_bg = int8(round( im_bg/num_frames_range ));


%% Analyse Images

%https://au.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html


for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    
    frame_input = im(:,:,frameNum);

    % ---- METHOD 1 ------------ %
%     %     frame_output = method_01(frame_input, PLOT_NUMS);
%     % Maximise contrast, edge detection, dilation and filling
%     
%     %Read image and modify the contrast
%     I1_orig = frame_input;
%     I1 = I1_orig;
%     I1 = histeq(I1_orig);
%     
%     %Edge detection
%     I2 = edge(I1,'canny',0.97); %max cont: 0.97
%     %     I2 = filledgegaps(I2, 10);
%     I2 = imfill(I2, 'holes');
%     %     I2 = bwareaopen(I2, 150);
%     
%     
%     % Idea: use fill edge gaps and Fill in initial canny output,
%     % if area of lines from canny / area of outline of fill < 1
%     % then fill failed and shape is not complete
%     % So redo canny with lower threshold until fill works
%     %
%     % But how can I limit this to just the motile cells? Do I need to?
%     
%     %Dilate image
%     dilationFactor = 5.5;
%     se90 = strel('line', dilationFactor, 90);
%     se0 = strel('line', dilationFactor, 0);
%     I3 = imdilate(I2, [se90 se0]);
%     
%     %Fill gaps
%     I3 = imfill(I3, 'holes');
%     I3 = bwareaopen(I3, 800);
%     
%     %Clear borders
%     %I5 = imclearborder(I4, 4);
%     
%     
%     %Smoothen
%     seD = strel('diamond',1);
%     I3 = imerode(I3, seD);
%     I3 = imerode(I3, seD);
%     cellOutline = bwperim(I3);
%     
%     %Extract centroids
%     regions = regionprops(I3);
%     regions_cell = struct2cell(regions);
%     cent_cell = regions_cell(2,:);
%     centroids = cell2mat(cent_cell');
%     
%     %Label connected regions
%     cc = bwconncomp(I3);
%     labeled = labelmatrix(cc);
%     I3 = label2rgb(labeled);
%     
%     %Overlay outline on original
%     I4 = I1_orig;
%     I4(cellOutline) = 255;
%     
%     % Maybe we need an adaptive dilation factor?
%     % Check if motile dicty filled in
%     % - Look for object with centroid nearest to motile centroid in
%     % previous frame.
%     % If not filled in, increase dilation a smidge and check again
    
    % --- METHOD 2 ------------ %
    
    % 1: Read image and normalise contrast
    I1_orig = frame_input;
    I1 = I1_orig;
    I1 = histeq(I1_orig);
    
    
    % 2: Extract brightest pixels
    I2 = I1;
    I2 = (I2>250); 
    
    % 3: Additional processing
    I3 = I2;
    % % Remove noise
    I3 = bwareaopen(I3, 75);
    % % Dilate image
    dilationFactor = 3; 
    SE = strel('disk',dilationFactor);
    I3 = imdilate(I3, SE);
    I3 = bwareaopen(I3, 500);  
    %Smoothen
%     seD = strel('diamond',1);
%     I3 = imerode(I3, seD);
%     I3 = imerode(I3, seD);
    cellOutline = bwperim(I3);
    cellSkeleton = bwperim(bwmorph(I3,'thin',Inf));
    
    
    %Extract centroids
    regions = regionprops(I3);
    regions_cell = struct2cell(regions);
    cent_cell = regions_cell(2,:);
    centroids = cell2mat(cent_cell');
    
    %Label connected regions
    cc = bwconncomp(I3);
    labeled = labelmatrix(cc);
    I3 = label2rgb(labeled);
    
    
    
    %Overlay outline on original
    I4 = I1_orig;
    I4(cellOutline) = 255;
    I4(cellSkeleton) = 255;
    
    
    % ------------
    
    frameIndex = frameNum - (FRAME_RANGE(1)-1);
    if frameNum == FRAME_RANGE(1)
        cent_hist = NaN * zeros(length(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)), 2);
        peak_hist = NaN * zeros(length(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)), 2);
        followIndex = 2; %HARDCODE follow initial specified centroid
        
    else
        cent_dists = sqrt(sum((centroids - cent_hist(frameIndex-1,:)).^2, 2));
        followIndex = find(cent_dists == min(cent_dists));
        
    end
    cent_hist(frameIndex,:) = centroids(followIndex, :);
    im_cell = (labeled == followIndex);
    
    [numPeaks, ~] = skelmetric(im_cell, 0.05);
    peak_hist(frameIndex,1) = numPeaks;
    [numPeaks, ~] = shapemetric(im_cell);
    peak_hist(frameIndex,2) = numPeaks;
    
    
    %Plot images
    plot_num_rows = PLOT_NUM_ROWS_LIST(num_plots);
    plot_num_cols = PLOT_NUM_COLS_LIST(num_plots);
    plot_num = 0;
    for im_num = PLOT_NUMS
        plot_num = plot_num + 1;
        plot_i = floor((plot_num - 1)/plot_num_cols);
        plot_j = mod(plot_num - 1, plot_num_cols);
        pos = zeros(1,4);
        pos(1) = plot_j/plot_num_cols + PLOT_MARGIN_X;
        pos(2) = 1 - (plot_i+1)/plot_num_rows + PLOT_MARGIN_Y;
        pos(3) = 1/plot_num_cols - 2*PLOT_MARGIN_X;
        pos(4) = 1/plot_num_rows - 2*PLOT_MARGIN_Y;
        subplot('Position', pos);
        eval(['imshow(I' num2str(im_num) ',''InitialMagnification'', ''fit'');']);
        title(IMAGE_NAMES(im_num), 'FontSize',16);
        
        if im_num == PLOT_NUMS(end)
            hold on;
            plot(centroids(:,1), centroids(:,2), 'rx');
            plot(cent_hist(:,1), cent_hist(:,2), 'b--');
            hold off;
            text(centroids(followIndex,1)+10, centroids(followIndex,2), ...
                [num2str(numPeaks) ' arms. '],...
                'Color', 'g');
        end
        
    end
    set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)], 'NumberTitle', 'off');
    drawnow();
    
    stepFrames = 0;
%     if (frameNum > 269) 
%         stepFrames = 1;
%         disp(['Pause on frame ' num2str(frameNum)]);
%     end
    if stepFrames == 1
        pause;
    end
    
end

figure(2);
x = 1:size(peak_hist,1);
y1 = peak_hist(:,1);
y2 = peak_hist(:,2);
plot(x, y1, 'b--x', x, y2, 'r--+');
legend('skelmetric', 'shapemetric');
title('comparison of both methods');

1;
end




% function frame_output = method_01_test(frame_input, PLOT_NUMS)
% Maximise contrast, edge detection, dilation and filling
%
%     Read image and modify the contrast
%     I1_orig = frame_input;
%     I1 = histeq(I1_orig);
%     I1 = (I1 > 247);
%
%     Edge detection
%     [~, threshold] = edge(I1,'sobel');
%     fudgeFactor = 2.28;
%     I2 = edge(I1, 'sobel', threshold * fudgeFactor);
%     I2 = edge(I1,'canny',0.925);
%
%     Dilate image
%     dilationFactor = 3.95;
%     se90 = strel('line', dilationFactor, 90);
%     se0 = strel('line', dilationFactor, 0);
%     I3 = imdilate(I2, [se90 se0]);
%
%     Fill gaps
%     I4 = imfill(I3, 'holes');
%     region = regionprops(I4);
%     Idea:
%     - loop over each connected region and see how much they differ from a
%     circle (center of bounding box, radius as midpoint of rectangle
%     edges).
%     - if the number of pixels that differ is less than a threshold, use
%     that disk to create a mask for the next frame
%
%     for i = 1:length(region)
%         cent = region(i).Centroid;
%         box = region(i).BoundingBox;
%         center = [box(1)+(box(3))/2, box(2)+(box(4))/2];
%     end
%
%
%     Clear borders
%     I5 = imclearborder(I4, 4);
%     I5 = bwareaopen(I4, 750);
%
%     Smoothen
%     seD = strel('diamond',1);
%     I6 = imerode(I5, seD);
%     I6 = imerode(I6, seD);
%
%     Overlay outline on original
%     cellOutline = bwperim(I6);
%     I7 = I1_orig;
%     I7(cellOutline) = 255;
%
%     frame_output = cell(1,max(PLOT_NUMS));
%     for i = PLOT_NUMS
%         eval(['frame_output{i} = I' num2str(i) ';'])
%     end
%
% end
%
%
%
% function [frame_output, data_output] = method_01(frame_input, PLOT_NUMS)
%
% end
%
%
% function frame_output = method_02_test(frame_input, PLOT_NUMS)
% Maximise contrast, threshold cell interiors, dilate
%
%     Maximise contrast and threshold to isolate cell interiors
%     I1_orig = frame_input;
%     I1 = histeq(I1_orig);
%     I1 = (I1 > 250);
%
%     Fill in image
%     I2 = imfill(I1, 'holes');
%
%     Remove Noise
%     I3 = bwareaopen(I2, 100);
%
%     Dilate image
%     I4 = imclose(I3, true(2));
%
%
%     %Smoothen
%     seD = strel('diamond',1);
%     I5 = imerode(I4, seD);
%     I5 = imerode(I5, seD);
%     I5 = I4;
%     I6 = I4;
%
%     Overlay outline on original
%     cellOutline = bwperim(I6);
%     I7 = I1_orig;
%     I7(cellOutline) = 255;
%
%     frame_output = cell(1,max(PLOT_NUMS));
%     for i = PLOT_NUMS
%         eval(['frame_output{i} = I' num2str(i) ';'])
%     end
%
% end
%
%
% function frame_output = method_02(frame_input, PLOT_NUMS)
%
%     Maximise contrast and threshold to isolate cell interiors
%     I1_orig = frame_input;
%     I1 = histeq(I1_orig);
%     I1 = (I1 > 238);
%
%     Fill in image
%     I2 = imfill(I1, 'holes');
%
%     Remove Noise
%     I2 = bwareaopen(I2, 250);
%
%     Dilate image
%     I2 = imclose(I2, true(7));
%
%     Overlay outline on original
%     cellOutline = bwperim(I2);
%     I3 = I1_orig;
%     I3(cellOutline) = 255;
%
%     frame_output = cell(1,max(PLOT_NUMS));
%     for i = PLOT_NUMS
%         eval(['frame_output{i} = I' num2str(i) ';'])
%     end
%
% end


% for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
%
%     Read image
%     I1 = im(:,:,frameNum);
%
%     Subtract background
%     I2 = int8(I1) - im_bg;
%     I2 = int8(I2>=0).*I2; %Interior
%     I2 = -int8(I2<=0).*I2; %Exterior
%     I2 = uint8(I2);
%
%     Test using threshold
%     I3 = zeros(ROI_height,ROI_width,'logical');
%     for i = 1:ROI_height
%         for j = 1:ROI_width
%             if I2(i,j) < 8
%                 I3(i,j) = 0;
%             else
%                 I3(i,j) = 1;
%             end
%         end
%     end
%
%     Remove small and large areas
%     I4 = I3 - bwareaopen(I3,AREA_LIMITS(2));
%     I4 = bwareaopen(I4,AREA_LIMITS(1));
%
%     %Edge detection
%     I4 = edge(I3, 'canny',0.37);
%
%     Try looking at the bright spots? background is grey, cell exterior is
%     black, cell interior is white..
%
%     Plot images
%     plot_num_rows = PLOT_NUM_ROWS_LIST(num_plots);
%     plot_num_cols = PLOT_NUM_COLS_LIST(num_plots);
%     plot_num = 0;
%     for im_num = PLOT_NUMS
%         plot_num = plot_num + 1;
%         plot_i = floor((plot_num - 1)/plot_num_cols);
%         plot_j = mod(plot_num - 1, plot_num_cols);
%         pos = zeros(1,4);
%         pos(1) = plot_j/plot_num_cols + PLOT_MARGIN_X;
%         pos(2) = 1 - (plot_i+1)/plot_num_rows + PLOT_MARGIN_Y;
%         pos(3) = 1/plot_num_cols - 2*PLOT_MARGIN_X;
%         pos(4) = 1/plot_num_rows - 2*PLOT_MARGIN_Y;
%         subplot('Position', pos);
%         eval(['imshow(I' num2str(im_num) ',''InitialMagnification'', ''fit'');']);
%         title(IMAGE_NAMES(im_num), 'FontSize',16);
%
%     end
%     set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)], 'NumberTitle', 'off');
%     drawnow();
% end



% %% Try edge dection?
% testframe = 200;
%
% im_sobel = edge(im(:,:,testframe), 'sobel');
% im_canny = edge(im(:,:,testframe), 'canny');
% im_prewitt = edge(im(:,:,testframe), 'prewitt');
%
% %% testplots
%
% figure;
% I1 = im(:,:,testframe);
% imshow(I1,'InitialMagnification','fit');
%
% figure;
% I2 = int16(im(:,:,testframe)) - im_bg(:,:);
% I2 = uint8(I2);
% imshow(I2,'InitialMagnification','fit');
%
% % figure;
% % I2 = im_canny(:,:);
% % for i = 1:size(REGIONS_TO_IGNORE,1)
% %     for j = REGIONS_TO_IGNORE(i,2):REGIONS_TO_IGNORE(i,4)
% %         for k = REGIONS_TO_IGNORE(i,1):REGIONS_TO_IGNORE(i,3)
% %             I2(j,k) = 0;
% %         end
% %     end
% % end
% % imshow(I2,'InitialMagnification','fit');




% IDEA: fill in areas close by centroids with low variance in position??
