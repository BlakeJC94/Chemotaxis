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


clipNum = 2;
makeMovie = 0;
frameRate = 10;
stepEvent = 0;


DIR = 'data/'; %B: Directory changed to work on my laptop
load_clip_and_pars;



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

if makeMovie == 1
    PLOT_NUMS = 1:2;
    IMAGE_NAMES = {'Original', 'Overlay'};
end


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


%% Analyse Images

if makeMovie == 1
    append = 0;
    movieName = ['Clip_' num2str(clipNum,'%02.f') '_' FILENAME(1:end-4) ...
        '_overlay_' num2str(append,'%03.f') '.avi'];
    
    while exist(['videos/' movieName], 'file') ~= 0 
        append = append+1;
        movieName = ['Clip_' num2str(clipNum,'%02.f') '_' FILENAME(1:end-4) ...
        '_overlay_' num2str(append,'%03.f') '.avi'];
    end
    
    v = VideoWriter(['videos/' movieName]); %create writer object
    v.FrameRate = frameRate;
    open(v); %open obj
end


watershedFlag = 0;

for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    
    frame_input = im(:,:,frameNum);

    
    % 1: Read image, remove non-uniform bg with morphopen and adjust contrast
    I1_orig = frame_input;
    I1 = I1_orig;
    background = imopen(I1, strel('disk',15));
    I1 = I1 - background;
    I1 = imadjust(I1);
    
    
    % 2: Get cell outline by threshold method and imclose
    I2 = im2bw(I1, graythresh(I1));
    I2 = imfill(I2,'holes');
    I2 = imclose(I2, ones(3)); %dilate then erode
    I2 = bwareaopen(I2, noiseThr); %remove false positives
    SE = strel('disk',3);
    I2 = imdilate(I2, SE); %dilate to enclose cells
    
    
    % 3: Seperate using internal markers watershed transform
    mask_em = imextendedmax(I1, 1); %find local maxima (excl. small max below)
    mask_em = bwareaopen(mask_em, noiseThr); %remove false positives and noise
    mask_em = imclose(mask_em, ones(4)); %dilate then erode
    mask_em = imfill(mask_em, 'holes');
    
    I_eq_c = imcomplement(I1);
    I_mod = imimposemin(I_eq_c, ~I2 | mask_em); %bg and interiors forced to be min
    L = watershed(I_mod);
    I3 = (L>1); %remove background and outlines
    I3 = bwareaopen(I3, noiseThr); %remove small noise
    cc = bwconncomp(I3);
    labeled = labelmatrix(cc);
    
    % 4: Overlay outline on original
    cellOutline = bwperim(I3);
    cellSkeleton = bwperim(bwmorph(I3,'thin',Inf));
    I4 = imoverlay(I1_orig, cellOutline | cellSkeleton, [.3 1 .3]);
    
    
    
    %Extract centroids
    regions = regionprops(I3);
    regions_cell = struct2cell(regions);
    cent_cell = regions_cell(2,:);
    centroids = cell2mat(cent_cell');
    
    %debug: Label connected regions from watershed
%     cc = bwconncomp(I3);
%     labeled = labelmatrix(cc);
    I3 = label2rgb(L);
 
    
    
    % ------------
    
    
    frameIndex = frameNum - (FRAME_RANGE(1)-1);
    if frameNum == FRAME_RANGE(1)
        cent_hist = NaN * zeros(length(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)), 2);
        peak_hist = NaN * zeros(length(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)), 2);
%         followIndex = 1; %HARDCODE follow initial specified centroid
        
    else
        cent_dists = sqrt(sum((centroids - cent_hist(frameIndex-1,:)).^2, 2));
        followIndex = find(cent_dists == min(cent_dists));
        if numel(cent_dists) > 1
            tmp = unique(cent_dists(:));
            next_cent_dist = tmp(2);
        end
        
    end
    cent_hist(frameIndex,:) = centroids(followIndex, :);
    im_cell = (labeled == followIndex);
    
    [numPeaks, ~] = skelmetric(im_cell, metricThreshold);
    peak_hist(frameIndex,1) = numPeaks;
    
    if makeMovie == 1
        I1 = I1_orig;
        I2 = I4;
    end
    
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
            text(10, 10, ...
                ['No. of arms : ' num2str(numPeaks)],...
                'Color', 'b');
        end
        
    end
    set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)], 'NumberTitle', 'off');
    drawnow();
    
    if makeMovie == 1
        frame = getframe(gcf); %save current figure as frame
        writeVideo(v,frame); %write frame to file
    end
    
    stepFrames = 0;
    if (stepEvent == 1) && (frameNum == 200)
            stepFrames = 1;
    end
    if stepFrames == 1
        str = ['Paused on frame ' num2str(frameNum) '. Press any key to continue \n'];
        fprintf(str);
%        1;
        pause;
        fprintf(repmat('\b',1,length(str)-1));
        stepFrames = 0;
    end
    
end

if makeMovie == 1
    close(v);
end


end







%     % 2: Extract brightest pixels
%     I2 = I1;
%     I2 = (I2>histeqThreshold); 
%     
%     
%     % 3: Additional processing
%     I3 = I2;
%     % % Remove noise
%     I3 = bwareaopen(I3, smallnoiseThreshold);
%     % % Dilate image
%     SE = strel('disk',dilationFactor);
%     I3 = imdilate(I3, SE);
%     % % Smoothen
%     seD = strel('diamond',1);
%     I3 = imerode(I3, seD);
%     I3 = imerode(I3, seD);
%     
%     % % Remove large outliers
%     I3 = bwareaopen(I3, largenoiseThreshold); 
%     I3 = imfill(I3,'holes'); 
%     
%     if watershedFlag == 1
%         
% %         % % --- method 1 ---- 
% %         % Identify foreground
% %         mask_em = (I1>seperationThreshold); 
% %         mask_em = bwareaopen(mask_em, 200);
% %         mask_em = imfill(mask_em,'holes'); 
% %         % Compute Watershed ridge lines and remove from I3 to seperate cells
% %         I3_D = bwdist(mask_em);
% %         I3_DL = watershed(I3_D);
% %         I3_bgm = (I3_DL == 0);
% %         I3(I3_bgm == 1) = 0;
%         
%         
%         % % --- method 2 ----
%         %Find internal markers (brighter spots)
%         mask_em = imextendedmax(I1, 1); %%
%         mask_em = bwareaopen(mask_em, 400); %%
%         mask_em = imclose(mask_em, ones(5,5));
%         mask_em = imfill(mask_em, 'holes');
%         
%         %complement image to make peaks into valleys
%         I_eq_c = imcomplement(I1);
%         
%         %modify image such that bg and maxima are the only local minima
%         I_mod = imimposemin(I_eq_c, ~I3 | mask_em);
%         
%         %compute watershed
%         L = watershed(I_mod);
%         I3_bgm = (L == 0);
%         I3_bgm = imdilate(I3_bgm,strel('disk',2));
%         I3(I3_bgm == 1) = 0;
%         
% %         % % ---- method 3 ----
% %         % Identify foreground
% %         mask_em = (I1>254); 
% %         mask_em = bwareaopen(mask_em, 300);
% %         mask_em = imfill(mask_em,'holes'); 
% %         %complement image to make peaks into valleys
% %         I_eq_c = imcomplement(I1);
% %         %modify image such that bg and maxima are the only local minima
% %         I_mod = imimposemin(I_eq_c, ~I3 | mask_em);
% %         %compute watershed
% %         L = watershed(I_mod);
% %         I3_bgm = (L == 0);
% %         I3(I3_bgm == 1) = 0;
% %         I3 = bwareaopen(I3, largenoiseThreshold); 
% %         I3 = imfill(I3,'holes'); 
%         
%     end
%     %Maybe watershed should only be active when two centroids are within
%     %range of overlapping?
    
