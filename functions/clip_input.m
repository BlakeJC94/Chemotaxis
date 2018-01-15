function [FILENAME, CHANNEL, FRAME_JUMP,FRAME_RANGE, ROI] = clip_input(clipNum)
%Stores arguemnts for AnalyseDictyShape, indexed by positive integers

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
    FRAME_RANGE = [-1 Inf];
    FRAME_JUMP = 1;
    ROI = [170, 200, 710, 500];
    
else
    error(['No clip corresponding to ' num2str(clipNum)]);
end

end

