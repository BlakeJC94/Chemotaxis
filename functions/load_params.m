function [histeqThreshold, smallnoiseThreshold, dilationFactor,largenoiseThreshold, followIndex, metricThreshold] = load_params(clipNum)
%Stores arguemnts for AnalyseDictyShape, indexed by positive integers
%    histeqThreshold: extract pixels with value > threshold (max 255)
%    smallnoiseThreshold: fills in holes with size less than the threshold
%    dilationFactor: replace all pixels with disks of radius dilationFactor
%    largenoiseThreshold: removes larger shapes cased by dilating noise
%    folowIndex: track centroid initially corresponding to followIndex

if clipNum == 1
    %Orignal test clip
    histeqThreshold = 180;
    smallnoiseThreshold = 20; 
    dilationFactor = 4; 
    largenoiseThreshold = 900; 
    followIndex = 2; 
    metricThreshold = 0.05;
    
elseif clipNum == 2
    histeqThreshold = 250;
    smallnoiseThreshold = 85; 
    dilationFactor = 5; 
    largenoiseThreshold = 900; 
    followIndex = 1; 
    metricThreshold = 0.05;
    
else
    error(['No clip corresponding to ' num2str(clipNum)]);
end

end

