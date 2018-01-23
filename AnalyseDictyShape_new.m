function AnalyseDictyShape_new
% Function for analysing the images from Dicty data
% Based on code from David M. Richards - 14/08/2017
% Blake Cook - 6/12/2017
% ======================== %

%% Preamble
%  - clipNum: integer index of clip to load (input for loading data script)
%  - makeMovie: (1) show original data with final overlay and export 
%      avi with specified |frameRate| (0) show full output of analysis and 
%      don't export avi
%  - stepEvent: [DEBUG] (1) pause analysis for some condition specified in
%      string |cond| (0) don't pause analysis until finished
% Load clip using script |load_clip_and_pars|

addpath('functions/') % add dir for functions
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning'); 
clf; % clear current figure
set(gcf, 'Position', get(0,'Screensize')); % set size of figure window

DIR = 'data/'; % specify where data is relative to current folder

makeMovie = 0;
frameRate = 10;
stepEvent = 0;

clipNum = 1;




%% Plot parameters:
%  - PLOT_NUMS: plot numbers to show
%  - PLOT_MARGIN: margins around each subplot
%  - PLOT_NUM_COLS_LIST, PLOT_NUM_ROWS_LIST: number of rows/rows for subplot
%  - IMAGE_NAMES: names for images

PLOT_NUMS = 1:6;
PLOT_MARGIN_X = 0.003;
PLOT_MARGIN_Y = 0;
PLOT_NUM_ROWS_LIST = [ 1 1 2 2 2 2 3 3 ];
PLOT_NUM_COLS_LIST = [ 1 2 2 2 3 3 3 3 ];
IMAGE_NAMES = {...
    'frame\_adjusted',...
    'ext\_binary',...
    'ext\_mask',...
    'int\_max',...
    'int\_mask',...
    'overlay'...
    };

% If makeMovie == 1, change plot output to original + overlay instead of
% full output from image analsysis. 
if makeMovie == 1
    PLOT_NUMS = 1:2;
    IMAGE_NAMES = {'Original', 'Overlay'};
end




%% Load clip and read images
% Use script to load clip and alansysis parameters stored in 
% |load_clip_andPars|: add the following variables to the workspace:
%     FILENAME, CHANNEL, FRAME_JUMP, FRAME_RANGE, ROI, noiseThr, metricThreshold
%  - im: 3d array that stores all images (rows, columns, frames)

% load_clip_and_pars; %TODO: FIX THIS!!!
FILENAME = 'DictyElectrotaxis_171116_001.tif';
CHANNEL = [3 3];
FRAME_JUMP = 1;
FRAME_RANGE = [155 271];
ROI = [340 500 880 800];

noiseThr = 100; %200

followIndex = 3;
metricThreshold = 0.05;



[im, num_plots, num_frames]  = readImages();


% If makeMovie == 1, specify avi output in 'Videos/' and set up video writer
if makeMovie == 1
    
    append = 0; % index to append to filename if duplicate is detected
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




%% Analyse images, get data and plot results
% to track cells, 

for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    
    frame_input = im(:,:,frameNum);
    analysisResults = analyseImages(frame_input);
    
    
    % get region data (centroids, area)
     %%%%%
    
    % update trackers
    if frameNum == FRAME_RANGE(1)
        
        % set up |cell_id|
        
    else 
        
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
        eval(['imshow(analysisResults.I' num2str(im_num) ...
            ', ''InitialMagnification'', ''fit'');']);
        title(IMAGE_NAMES(im_num), 'FontSize',16);
        
%         if im_num == 6
%             hold on;
% %             plot(centroids(:,1), centroids(:,2), 'rx');
% %             x3 = cent_hist{3}(:,2);
% %             y3 = cent_hist{3}(:,3);
% %             plot(x3,y3,'b--');
%             %             plot(cent_hist(:,1), cent_hist(:,2), 'b--');
%             hold off;
%             text(10, 10, ...
%                 ['No. of arms : ' num2str(numPeaks)],...
%                 'Color', 'b');
%             
%             for j = 1:size(active_id,2)
%                 text(centroids(active_id(2,j),1)-10,...
%                     centroids(active_id(2,j),2)-35, ...
%                     num2str(active_id(1,j)), 'FontSize', 18);
%             end
%             
%         end
        
    end
    set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)], 'NumberTitle', 'off');
    drawnow();
    
    
end



















%% SUBFUNCTION 1: readImages()
% Takes plot options specified in |AnalyseDictyShape| and analysis
% parameters stored in |load_clip_and_pars|.
% Required variables: DIR, FILENAME, CHANNEL, FRAME_RANGE, FRAME_JUMP, ROI,
%    PLOT_NUMS.

    function [im, num_plots, num_frames] = readImages()
        % Read header info
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
        
        % Declare variables:
        %  - im: holds all images
        im = zeros(ROI_height,ROI_width,num_frames,'uint8');
        
        % Read in images
        tif_obj = Tiff([DIR FILENAME],'r');
        for i = 1:num_frames
            tif_obj.setDirectory(CHANNEL(1)+CHANNEL(2)*(i-1));
            im_full = tif_obj.read();
            im(:,:,i) = im_full(ROI(2):ROI(4),ROI(1):ROI(3));
        end
        tif_obj.close();
        
    end




%% SUBFUNCTION 2: analyseImages(frame_input)
% Output: 
%  - I1: Adjust frame noise and contrast levels 
%  - I2: Threshold I1 for finding external mask
%  - I3: Get external mask of cells from I2 (do not worry about overlapping)
%  - I4: Find local maxima for finding internal mask
%  - I5: Get internal mask of cells from I4
%  - I6: Compute watershed transform and overlay on original image
% 
% The technique I utilise here uses the watershed transform. The MATLAB
% documentation explains this transform in the following way: imagine that 
% bright areas are "high" and dark areas are "low" and flow the flow lines 
% that water would naturally follow on this surface. regions can be
% seperated by using the "catchement basins", local minumms on the surface.
% 
% Our strategy involves extracting internal markers from the frame and
% forcing these to be catchement basis. 

    function analysisResults = analyseImages(frame_input)
        
        analysisResults = struct;
        
        % remove non-uniform noise from background
        bg = imopen(frame_input, strel('disk',15));
        noise_correct = frame_input - bg;
        % adjust constrast (1% of data saturated at low and high intensities)
        frame_adjusted = imadjust(noise_correct);
        
        I1 = frame_adjusted;
        
        
        % convert to binary image using threshold calculated to seperate
        % black and white pixels via |graythresh|
        ext_binary = imbinarize(frame_adjusted, graythresh(frame_adjusted));
        
        I2 = ext_binary;
        
        
        % fill holes (areas surronded by white pixels)
        ext_fill = imfill(ext_binary, 'holes');
        
        % morphologially close image (dilate then erode, closes gaps
        % between some regions) and remove false positives
        ext_fill_close = imclose(ext_fill, ones(3)); 
        ext_fill_close = bwareaopen(ext_fill_close, noiseThr); %TODO: define new var for this thr
        
        % dialate regions by disk (r=3 pix) to enclose cells
        ext_mask = imdilate(ext_fill_close, strel('disk',3));
        
        I3 = ext_mask;
        
        
        
        % find local maxima in image
        int_max = imextendedmax(frame_adjusted, 1);
        
        I4 = int_max;
        
        
        
        % remove noise
        int_max_main =  bwareaopen(int_max, noiseThr); 
        
        % morphologically close (dialte then erode, close some gaps) then
        % fill holes
        int_close = imclose(int_max_main, ones(4));
        int_mask = imfill(int_close, 'holes');
        
        I5 = int_mask;
        
        
        
        % complement |frame_adjust|
        im_c = imcomplement(frame_adjusted);
        
        % modify image such that background and internal mask are minumum
        % points in the image
        im_mod = imimposemin(im_c, ~ext_mask | int_mask);
        
        % computer watershed transform on |im_mod| and remove bg and
        % outlines of cells, leving distinct regions for cells
        L = watershed(im_mod);
        L = (L>1);
        
        % remove any noise introduced by watershed (outlines sometimes
        % fragment very small parts of the cell membrane)
        mask = bwareaopen(L, noiseThr);
        
        I6 = mask;
        
        
        
        % overlay cell boundaries and skeletons on original frame
        cellOutline =  bwperim(mask);
        cellSkeleton = bwperim(bwmorph(mask,'thin',Inf));
        overlay = imoverlay(frame_input, cellOutline | cellSkeleton, [.3 1 .3]);
        
        analysisResults.overlay = overlay;
        
        
        
        if makeMovie == 1
            I1 = frame_input;
            I2 = I6;
        end
        
        % export I1, I2,... to function output struct
        for i = PLOT_NUMS
            eval(['analysisResults.I' num2str(i) ' = I' num2str(i) ';'])
        end
    end

end




