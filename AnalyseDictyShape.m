function AnalyseDictyShape
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
%      string |stepCond| (0) don't pause analysis until finished
% Load clip using script |load_clip_and_pars|

addpath('functions/') % add dir for functions
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning'); 
clf; % clear current figure
set(gcf, 'Position', get(0,'Screensize')); % set size of figure window

DIR = 'data/'; % specify where data is relative to current folder

makeMovie = 0;
frameRate = 10;

stepEvent = 0;
stepCond = 'frameNum == 200';

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
% Use script to load clip and analysis parameters stored in 
% |load_clip_andPars|: add the following variables to the workspace:
%     FILENAME, CHANNEL, FRAME_JUMP, FRAME_RANGE, ROI, noiseThr, metricThreshold
%  - im: 3d array that stores all images (rows, columns, frames)

% initialise variables
FILENAME = []; CHANNEL = []; FRAME_JUMP = []; FRAME_RANGE = []; ROI = [];
noiseThr = []; followIndex = []; metricThreshold = [];

% load values
load_clip_and_pars; 

% read images 
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
% to track cells, ....
% Centroids matrix is indexed by a region index assigned by MATLAB. "Update
% trackers" segment matches these region indicies with the cell ID. 
for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    
    frame_input = im(:,:,frameNum);
    analysisResults = analyseImages(frame_input);
    
    
    % get region data (area, centroid) 
    regions = regionprops(analysisResults.I7, {'Centroid', 'Area'});
    regions = struct2cell(regions);
    
    areas = cell2mat(regions(1,:)');
    centroids = cell2mat(regions(2,:)');
    
    % get armnum data
    armnums = skelmetric_mult(analysisResults.I7, metricThreshold);
    
    % update ID and trackers
    if frameNum == FRAME_RANGE(1)
        
        % set up |cell_id|. 1st row = cell ID, 2nd row = region index (RI) 
        % in current frame.
        id = [(1:size(centroids,1)); (1:size(centroids,1))];
        ind = (id(2,:)~=-1);
        active_id = id(:,ind);
        
        % set up centroid tracker: cent_hist{n} returns a matrix of the
        %   centroid loactions of cell with ID 'n' (columns are frame, x, y)
        % set up area tracker
        % set up armnum trackercl
        cent_hist = cell(1,size(id,2));
        area_hist = cell(1,size(id,2));
        arm_hist = cell(1,size(id,2));
        for j = 1:size(id,2)
            cent_hist{id(1,j)} = [1, centroids(id(2,j),:)];
            area_hist{id(1,j)} = [1, areas(id(2,j))];
            arm_hist{id(1,j)} = [1, armnums(id(2,j))];
        end
        
        
        
    else 
        
        % update region index for each cell ID
        %  - assignment: nx2 matrix, col 1 = RIs on previous frame, 
        %     col 2 = corresponding RIs on current frame.
        %  - unassignedTracks: previous RIs not assigned (exited cells)
        %  - unassignedDetections: current RIs new (entering cells)
        cost = zeros(size(centroids_old,1), size(centroids,1));
        for j = 1:size(centroids_old,1)
            diff = centroids - repmat(centroids_old(j,:), [size(centroids,1),1]);
            cost(j,:) = sqrt(sum(diff.^2,2));
        end
        [assignment,unassignedTracks,unassignedDetections] = ...
            assignDetectionsToTracks(cost, 50);
        
        % create new id matrix
        id_update = [id(1,:); zeros(1,length(id(2,:)))];
        
        % disable IDs for cells that exit
        id_update(id == -1) = -1; % for cells that have already exited
        if ~isempty(unassignedTracks)
            disp(['unassignedTracks, frameNum = ' num2str(frameNum)])
            ind = ismember(id(2,:), unassignedTracks(:));
            id_update(2,ind) = -1;
        end
        
        
        % update RIs (assignment: col 1 prev, col 2 current)
        for l = 1:size(assignment,1)
            ind = ismember(id(2,:),assignment(l,1));
            id_update(2,ind) = assignment(l,2);
        end
        
        % add IDs for cells that enter
        if ~isempty(unassignedDetections)
            disp(['unassignedDetections, frameNum = ' num2str(frameNum)])
            
            % find last recorded positions of inactive cells
            inactive_id = id(1,(id(2,:) == -1));
            inactive_cents = zeros(length(inactive_id), 2);
            for l = 1:length(inactive_id)
                inactive_cents(l,:) = cent_hist{inactive_id(l)}(end,2:3);
            end
            
            
                
            for k = 1:length(unassignedDetections)
                
                % check if this cell already has an ID
                cent_dist = sqrt(sum((centroids(unassignedDetections(k),:) - inactive_cents).^2,2));
                if min(cent_dist) <= 20
                    cell_id = inactive_id((cent_dist == min(cent_dist)));
                    id_update(2, id(1,:) == cell_id) = unassignedDetections(k);
                else
                    new_id = max(id_update(1,:)) + 1;
                    id_update = [id_update(1,:), new_id; id_update(2,:), double(unassignedDetections(k))];
                
                    cent_hist{new_id} = zeros(0); 
                    area_hist{new_id} = zeros(0); 
                    arm_hist{new_id} = zeros(0);
                end
            end
        end
        
        % update ID matrix 
        id = id_update;
        
        % update trackers
        ind = (id(2,:)~=-1);
        active_id = id(:,ind);
        for j = 1:size(active_id,2)
            cent_hist{active_id(1,j)} = ...
                [cent_hist{active_id(1,j)}; frameNum, centroids(active_id(2,j),:)];
            
            area_hist{active_id(1,j)} = ...
                [area_hist{active_id(1,j)}; frameNum, areas(active_id(2,j))];
            
            arm_hist{active_id(1,j)} = ...
                [arm_hist{active_id(1,j)}; frameNum, armnums(active_id(2,j))];
        end
        
        
    end
    
    % store centroids for next frame 
    centroids_old = centroids;
    
    
    
    
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
        
        if im_num == PLOT_NUMS(end)
%             imshow(analysisResults.I7, 'InitialMagnification', 'fit');
            hold on;
            plot(centroids(:,1), centroids(:,2), 'rx');
            xn = cent_hist{followIndex}(:,2);
            yn = cent_hist{followIndex}(:,3);
            plot(xn,yn,'b--');
            hold off;
            text(10, 10, ...
                ['No. of arms : ' num2str(arm_hist{followIndex}(end,2))],...
                'Color', 'b');
            
            for j = 1:size(active_id,2)
                text(centroids(active_id(2,j),1)-10,...
                    centroids(active_id(2,j),2)-35, ...
                    num2str(active_id(1,j)), 'FontSize', 18);
            end
            
        end
        
    end
    % set title of figure and print results to figure now
    set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)], 'NumberTitle', 'off');
    drawnow();
    
    
    % if makeMovie == 1, export data to file. DEBUG: pause if 
    % stepEvent == 1 and not making movie
    if makeMovie == 1
        frame = getframe(gcf); %save current figure as frame
        writeVideo(v,frame); %write frame to file
        
    elseif (stepEvent == 1) && eval(stepCond)
        str = ['Paused on frame ' num2str(frameNum) '. Press any key to continue \n'];
        fprintf(str);
        pause; fprintf(repmat('\b',1,length(str)-1));
        
    end
    
    
end

% close video writer if opened
if makeMovie == 1
    close(v);
end


%% Plot results (RMS displacement)

% !! CLEAN THIS UP!



dataPoints = floor(numel(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2))/2);
plotData = cell(1,size(cent_hist,2));
index = 0;
for i = 1:size(cent_hist,2)
    if size(cent_hist{i},1) >= dataPoints
        index = index+1;
        %get dispalcement
        plotData{index} = sqrt(sum((cent_hist{i}(1:dataPoints,2:3) - cent_hist{i}(1,2:3)).^2,2));
    end
end

figure;
hold on; 
avgDisplacement = zeros(dataPoints,1);
for i = 1:index 
%     plot(1:dataPoints, plotData{i}(:), 'b--');
    avgDisplacement = ((i-1)*avgDisplacement + plotData{i}(:))/i;
end
plot(0:dataPoints-1, sqrt(avgDisplacement(:)), 'r--')
hold off;
title(['RMS displacement, cells present for at least ' num2str(dataPoints) ' frames']);



    

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
% forcing these to be catchement basis. Ueful figures and further
% explanation of watershed: 
% https://uk.mathworks.com/company/newsletters/articles/the-watershed-transform-strategies-for-image-segmentation.html
% 
% 
% Issues unresolved: sometimes internal mask will randomly fragment,
% resulting in incorrectly identifying one cell as 2 seperate cells 
%   - Possible fix: identify corresponding regions between internal mask
%   and external mask. If one region in external mask has more than one
%   region in internal mask, check if the area of exeternal mask has change
%   significantly. If so, overlapping has occurred and internal region
%   should be kept. Otherwise just take the largest internal marker. 


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
        
        % morphologially opemn image (erode then dilate, opens gaps
        % between some regions) and remove false positives
        ext_fill_close = imopen(ext_fill, ones(3)); 
        ext_fill_close = bwareaopen(ext_fill_close, noiseThr); %TODO: define new var for this thr
        
        % dialate regions by disk (r=3 pix) to enclose cells
        ext_mask = imdilate(ext_fill_close, strel('disk',3));
        
        I3 = ext_mask;
        
        
        
        % find local maxima in image
        int_max = imextendedmax(frame_adjusted, 1);
        
        I4 = int_max;
        
        
        
        % remove noise
        int_max_main =  bwareaopen(int_max, 75); 
        
        % morphologically close (dialte then erode, close some gaps), fill 
        % holes 
        int_close = imclose(int_max_main, ones(3));
        int_close_fill = imfill(int_close, 'holes');
%         int_thin = bwmorph(int_close_fill, 'Thin', 1);
        
        
        % remove false positives
        int_mask =  bwareaopen(int_close_fill, 25); 
        
        %%%% test defragmentation.. couldnt get this to work, may come back
        %%%% to this someday..
%         int_cent = cell2mat(struct2cell(regionprops(int_mask, 'Centroid'))');
%         int_cc = bwconncomp(int_mask);
%         int_labeled = labelmatrix(int_cc);
%         for i = 1:size(int_cent,1)
%             int_dists = sqrt(sum((int_cent(i,:)-int_cent).^2,2));
%             reg_del = find((int_dists < 150)&(int_dists > 0));
%             if ~isempty(reg_del)
%                 int_labeled(ismember(int_labeled,reg_del)) = 0;
%             end
%         end
%         
        
        
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
        
        I7 = mask;
        
        
        
        % overlay cell boundaries and skeletons on original frame
        cellOutline =  bwperim(mask);
        cellSkeleton = bwperim(bwmorph(mask,'thin',Inf));
        overlay = imoverlay(frame_input, cellOutline | cellSkeleton, [.3 1 .3]);
        
        I6 = overlay;
        
        
        
        if makeMovie == 1
            I1 = frame_input;
            I2 = I6;
        end
        
        % export I1, I2,... to function output struct
        for i = 1:7
            eval(['analysisResults.I' num2str(i) ' = I' num2str(i) ';'])
        end
        
    end


end



% function [armnums] = skelmetric_mult_old(im, peakThreshold)
% % Performs skelmetric transform on binary image of multiple regions. 
% %THIS IS AN INEFFICIENT WAY OF DOING THINGS
% cc = bwconncomp(im);
% labeled = labelmatrix(cc);
% 
% armnums = zeros(cc.NumObjects,1);
% for regionIndex = 1:cc.NumObjects
%     im_cell = (labeled == regionIndex); 
%     armnums(regionIndex) = skelmetric(im_cell, peakThreshold);
% end
% 
% end


function armnums = skelmetric_mult(im, peakThreshold)
% Performs skelmetric transform on binary image of multiple regions



% find skeletons
im_sm = bwmorph(im,'thin',Inf);
   
    
% find and dialte branchpoints
im_br = bwmorph(im_sm,'branchpoints');
dilationFactor = 2;
SE = strel('disk',dilationFactor);
im_br = imdilate(im_br, SE);
% remove from skeleton to get all skeletal connections
im_skelconn = logical(im_sm .* (1-im_br));

    
% find endpoints of skeleton
im_en = bwmorph(im_sm,'endpoints');
    
% invert nodeless skel, get br-br conn..
im_tmp = logical(1-im_skelconn);
% fill in br-en connections and invert to get br-br connections
[i,j] = find(im_en);
im_tmp = imfill(im_tmp, [i j], 8);
im_brconn = logical(1-im_tmp);
% subtract from all connections to obtain br-en connections
im_enconn = logical(im_skelconn - im_brconn);
    

cc = bwconncomp(im_sm);
labeled = labelmatrix(cc);
skelLengths = cell2mat(struct2cell(regionprops(im_sm,'Area')));

armnums = zeros(cc.NumObjects,1);
for currentcell = 1:cc.NumObjects
    
    % select skeleton and br-en connections
    curr_skeleton = im_sm(labeled == currentcell);
    curr_enconn = im_enconn(labeled == currentcell);
    
    
    curr_skelLength = skelLengths(currentcell);
    
    regions_enconn = regionprops(curr_enconn);
    armNum = 0;
    curr_skelLength = skelLengths(currentcell);
    if length(regions_enconn) == 1
        armNum = 2;
    else
        for k = 1:length(regions_enconn)
            area_en = regions_enconn(k).Area;
            if area_en/curr_skelLength > peakThreshold
                armNum = armNum + 1;
            end
        end
    end
    % if shape is a line, count as two arms
    if armNum == 1, armNum = 2; end
    
    armnums(currentcell) = armNum;
    
end



end






