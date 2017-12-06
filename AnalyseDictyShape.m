function AnalyseDictyShape
% Function for analysing the images from Dicty data
% Based on code from David M. Richards - 14/08/2017
% Blake Cook - 6/12/2017
% ======================== %

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
clf;
set(gcf, 'Position', get(0,'Screensize'));

% *** File "DictyElectrotaxis_171116_001" - 953Mb - 271 frames - 3 channels - 1280x960
DIR = 'data/'; %B: Directory changed to work on my laptop
FILENAME = 'DictyElectrotaxis_171116_001.tif';
CHANNEL = [3 3];
FRAME_RANGE = [155 157]; %B: looking at few frames for now, default "[155 271]"
FRAME_JUMP = 1;
ROI = [340 500 880 800];
REGIONS_TO_IGNORE = [95 5 185 100; 190 1 280 25; 320 1 350 20; 310 140 390 210; 480 225 510 260];



end