% ============================== %
% Analyse Dicty shape and motion %
% David M. Richards - 14/08/2017 %
% ============================== %

%Analysis parameters:
%  - CHANNEL - use channel number CHANNEL(1) out of total of CHANNEL(2) channels
%  - FRAME_RANGE - which frames to analysis; [-1 Inf] => all frames
%  - ROI - region of interest as [x1 y1 x2 y2]; [-1 -1 Inf Inf] => full image
%  - REGIONS_TO_IGNORE - regions to remove from final image as [x1 y1 x2 y2]

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
clearvars -except -regexp save$;
clf;
set(gcf, 'Position', get(0,'Screensize'));

% ------------------- %
% Analysis parameters %
% ------------------- %

% *** File "DictyElectrotaxis_171116_001" - 953Mb - 271 frames - 3 channels - 1280x960
DIR = 'data/'; %B: Directory changed to work on my laptop
FILENAME = 'DictyElectrotaxis_171116_001.tif';
CHANNEL = [3 3];
FRAME_RANGE = [155 185]; %B: looking at few frames for now, default "[155 271]"
FRAME_JUMP = 1;
ROI = [340 500 880 800];
REGIONS_TO_IGNORE = [95 5 185 100; 190 1 280 25; 320 1 350 20; 310 140 390 210; 480 225 510 260];

THRESHOLD = 8;
AREA_LIMITS = [500 10000];

% --------------- %
% Plot parameters %
% --------------- %

PLOT_NUMS = 1:7; %plot numbers to show
PLOT_MARGIN_X = 0.003; %x-margin to put around each subplot
PLOT_MARGIN_Y = 0; %y-margin to put around each subplot
PLOT_NUM_COLS_LIST = [ 1 2 2 2 3 3 3 3 ]; %number of columns for subplot
PLOT_NUM_ROWS_LIST = [ 1 1 2 2 2 2 3 3 ]; %number of rows for subplot
IMAGE_NAMES = {'Original image', 'Background subtracted', 'Thresholded', 'Holes filled', 'Small and large areas removed', 'Vanishing difference to background', 'Other cells ignored'}; %names for images

% ---------------- %
% Read header info %
% ---------------- %

file_info = imfinfo([DIR FILENAME]);
im_width = file_info(1).Width;
im_height = file_info(1).Height;
num_images = numel(file_info);

if mod(num_images,CHANNEL(2))~=0
    error(['Total number of frames (' num2str(num_images) ') not multiple of number of channels (' num2str(CHANNEL(2)) ')']);
end
num_frames = num_images/CHANNEL(2);

if (FRAME_RANGE(1)==-1); FRAME_RANGE(1)=1; end;
if (FRAME_RANGE(2)==Inf); FRAME_RANGE(2)=num_frames; end;

num_frames_range = numel(FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2));

if (ROI(1)==-1); ROI(1)=1; end;
if (ROI(2)==-1); ROI(2)=1; end;
if (ROI(3)==Inf); ROI(3)=im_width; end;
if (ROI(4)==Inf); ROI(4)=im_height; end;

ROI_width  = ROI(3) - ROI(1) + 1;
ROI_height = ROI(4) - ROI(2) + 1;

num_plots = numel(PLOT_NUMS);

% ----------------- %
% Declare variables %
% ----------------- %

im = zeros(ROI_height,ROI_width,num_frames,'uint8'); %holds all images
im_full = zeros(im_height,im_width); %holds full frame

% -------------- %
% Read in images %
% -------------- %

tif_obj = Tiff([DIR FILENAME],'r');
for i = 1:num_frames
   tif_obj.setDirectory(CHANNEL(1)+CHANNEL(2)*(i-1));
   im_full = tif_obj.read();
   im(:,:,i) = im_full(ROI(2):ROI(4),ROI(1):ROI(3));
end
tif_obj.close();

% --------------------- %
% Find background image %
% --------------------- %

im_bg = zeros(ROI_height,ROI_width,'double');
for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    for i = 1:ROI_height
        for j = 1:ROI_width
            im_bg(i,j) = im_bg(i,j) + double(im(i,j,frameNum));
        end
    end
end
im_bg = int16(round( im_bg/num_frames_range ));

% -------------- %
% Analyse images %
% -------------- %

for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)

    %Read image
    I1 = im(:,:,frameNum);

    %Subtract background
    I2 = abs( int16(I1) - im_bg );
    I2 = uint8(I2);

    %Threshold
    I3 = zeros(ROI_height,ROI_width,'logical');
    for i = 1:ROI_height
        for j = 1:ROI_width
            if I2(i,j) < THRESHOLD
                I3(i,j) = 0;
            else
                I3(i,j) = 1;
            end
        end
    end

    %Fill holes
    I4 = imfill(I3,'holes');

    %Remove small and large areas
    I5 = I4 - bwareaopen(I4,AREA_LIMITS(2));
    I5 = bwareaopen(I5,AREA_LIMITS(1));

    %Find white regions where difference to background vanishes
    I6 = zeros(ROI_height,ROI_width,'logical');
    for i = 1:ROI_height
        for j = 1:ROI_width
            if I5(i,j)==1 && I3(i,j)==0
                I6(i,j) = 1;
            end
        end
    end

    %Ignore some regions
    I7 = I6;
    for i = 1:size(REGIONS_TO_IGNORE,1)
        for j = REGIONS_TO_IGNORE(i,2):REGIONS_TO_IGNORE(i,4)
            for k = REGIONS_TO_IGNORE(i,1):REGIONS_TO_IGNORE(i,3)
                I7(j,k) = 0;
            end
        end
    end

    %Plot images
    plot_num_cols = PLOT_NUM_COLS_LIST(num_plots);
    plot_num_rows = PLOT_NUM_ROWS_LIST(num_plots);
    plot_num = 0;
    for im_num = PLOT_NUMS
        plot_num = plot_num + 1;
        plot_i = floor((plot_num-1)/plot_num_cols); plot_j = mod(plot_num-1, plot_num_cols);
        pos = [plot_j/plot_num_cols 1-(plot_i+1)/plot_num_rows 1/plot_num_cols 1/plot_num_rows] + [PLOT_MARGIN_X PLOT_MARGIN_Y -2*PLOT_MARGIN_X -2*PLOT_MARGIN_Y];
        subplot('Position',pos);

        eval(['imshow(I' num2str(im_num) ',''InitialMagnification'',''fit'');']);
        title(IMAGE_NAMES(im_num),'FontSize',16);
    end
    set(gcf,'Name',[FILENAME ': frame ' num2str(frameNum) '/' num2str(num_frames)],'NumberTitle','off');
    drawnow();
end