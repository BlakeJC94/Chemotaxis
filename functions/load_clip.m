function [FILENAME, CHANNEL, FRAME_JUMP,FRAME_RANGE, ROI] = load_clip(clipNum)
%Stores arguments for AnalyseDictyShape, indexed by positive integers

if clipNum == 1
    %Orignal test clip
    FILENAME = 'DictyElectrotaxis_171116_001.tif';
    CHANNEL = [3 3];
    FRAME_JUMP = 1;
    FRAME_RANGE = [155 271]; 
    ROI = [340 500 880 800];
    
elseif clipNum == 2
    FILENAME = 'DictyElectrotaxis_171116_002.tif';
    CHANNEL = [3 3];
    FRAME_RANGE = [120 Inf];
    FRAME_JUMP = 1;
    ROI = [170, 200, 710, 500];
        
elseif clipNum == 3
    %seems to be in negative, will require imcomplement() before analysis
    FILENAME = 'DictyElectrotaxis_171108_001.tif';
    CHANNEL = [1 3];
    FRAME_RANGE = [90 165]; %66 165
    FRAME_JUMP = 1;
    ROI = [170  362  710  662];
    
    
else
    error(['No clip corresponding to ' num2str(clipNum)]);
end

end

