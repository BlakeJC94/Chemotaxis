function ViewDictyData
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
% FILENAME = 'DictyElectrotaxis_171116_001.tif';
FILENAME = 'electrotaxis.tif';
CHANNEL = [2 3];
% FRAME_RANGE = [155 271]; 
FRAME_RANGE = [-1 Inf]; 
FRAME_JUMP = 1;
% ROI = [340 500 880 800];
ROI = [-1 -1 Inf Inf]; %Format: x1 y1 x2 y2

Rect = [10, 300, 800, 400]; %Format: x1 y1 w h 
ROC = [Rect(1), Rect(2), Rect(1)+Rect(3), Rect(2)+Rect(4)];


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

%https://au.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html


for frameNum = FRAME_RANGE(1):FRAME_JUMP:FRAME_RANGE(2)
    
    I1 = im(:,:,frameNum);
%     I1 = histeq(imcomplement(I1));
%     I1 = histeq(I1);
    imshow(I1,'InitialMagnification', 'fit');
    hold on;
    rectangle('Position',Rect,...
         'LineWidth',2,'LineStyle','--','EdgeColor','g');
    hold off; 
    
    title(['ROC : [' num2str(ROC) ']']);
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





end