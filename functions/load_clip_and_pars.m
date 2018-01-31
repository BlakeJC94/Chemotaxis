%Analysis parameters:
%  - CHANNEL - use channel number CHANNEL(1) out of total of CHANNEL(2) channels
%  - FRAME_RANGE - which frames to analysis; [-1 Inf] => all frames
%  - ROI - region of interest as [x1 y1 x2 y2]; [-1 -1 Inf Inf] => full image
%  - REGIONS_TO_IGNORE - regions to remove from final image as [x1 y1 x2 y2]
%  - noiseThr: 

if clipNum == 1
    %Orignal test clip
    FILENAME = 'DictyElectrotaxis_171116_001.tif';
    CHANNEL = [3 3];
    FRAME_JUMP = 1;
%     FRAME_RANGE = [155 271]; 
    FRAME_RANGE = [155 271];
%     ROI = [340 500 880 800];
    ROI = [-1 -1 Inf Inf];

    noiseThr = 100; %200
    
    followIndex = 2;  
    metricThreshold = 0.05;
    
elseif clipNum == 2
    FILENAME = 'DictyElectrotaxis_171116_002.tif';
    CHANNEL = [3 3];
    FRAME_RANGE = [120 Inf];
    FRAME_JUMP = 1;
    ROI = [170, 200, 710, 500];
    
    noiseThr = 10; 
    
    followIndex = 1;  
    metricThreshold = 0.05;
        
elseif clipNum == 3
    %seems to be in negative, will require imcomplement() before analysis
    FILENAME = 'DictyElectrotaxis_171108_001.tif';
    CHANNEL = [1 3];
    FRAME_RANGE = [90 165]; %66 165
    FRAME_JUMP = 1;
    ROI = [170  362  710  662];
    
    noiseThr = 200; 
    
    followIndex = 1;  
    metricThreshold = 0.05;
    
else
    error(['No clip corresponding to ' num2str(clipNum)]);
end