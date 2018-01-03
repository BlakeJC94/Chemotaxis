function [im_enconn, peakNum] = skelmetric(im)
% Performs skelmetric transform on binary image of single region. 
%
% Region is skeletonized and branchpoints and endpoints are identified.
% Branchpoints are dilated and removed from skeleton, leaving connections.
% branch-branch connections are removed, leaving only branch-endpoint
% connections.
% 
% Resulting image is used to determine how many 'arms' an irregular polygon
% has, which are counted as long branch-endpoint connections (relative to
% the total skeleton length), given as peakNum.
%
% Function will accept binary image of single region and return skelmetric
% binary image of br-en connections, as well as how many 'arms' the polygon
% has. 
%
% INPUTS:
%   im - binary image of a region
%
% OUTPUTS:
%   im_enconn - binary image of skeletal br-en connections
%   peakNum - number of 'arms' the polygon has
%
%
% Blake Cook 2017-12-27
% ------------------------

%get region data
regions = regionprops(im);
%if there is more than one region, select region with largest area
if length(regions)>1
    disp('More than one region detected, selection largest connected component');
    im_label = labelmatrix(bwconncomp(im));
    
    regions_cell = struct2cell(regions);
    areas = cell2mat(regions_cell(1,:));
    area_index = find(areas == max(areas));
    
    im = ismember(im_label, area_index) > 0;
    regions = regionprops(im);
end

%determine how round the object is
shapeArea = regions.Area;
regions_boundary = regionprops(bwperim(im));
shapePerimeter = regions_boundary.Area;
metric = 4*pi*shapeArea/shapePerimeter^2; %=1 iff shape is a circle 

%if the object is round, return msg and original image
%otherwise proceed with skeletonization
roundThreshold = 0.9;
if metric > roundThreshold
    
    disp(['shape is approximately round, metric = ' num2str(metric)]);
    im_enconn = im;
    peakNum = 0;
    
else
    
    %find skeleton
    im_sm = bwmorph(im,'thin',Inf);
    
%     %--debug--
%     subplot(2,2,2)
%     imshow(im_sm, [], 'InitialMagnification', 'fit');
%     title('skel');
    
    %find and dialte branchpoints
    im_br = bwmorph(im_sm,'branchpoints');
    dilationFactor = 2;
    SE = strel('disk',dilationFactor);
    im_br = imdilate(im_br, SE);
    %remove from skeleton to get all skeletal connections
    im_skelconn = logical(im_sm .* (1-im_br));
    
%     %--debug--
%     subplot(2,2,3)
%     imshow(im_skelconn, [], 'InitialMagnification', 'fit');
%     title('skel no nodes');
    
    %find endpoints of skeleton
    im_en = bwmorph(im_sm,'endpoints');
    
    %invert nodeless skel, get br-br conn..
    im_tmp = logical(1-im_skelconn);
    %fill in br-en connections and invert to get br-br connections
    [i,j] = find(im_en);
    im_tmp = imfill(im_tmp, [i j], 8);
    im_brconn = logical(1-im_tmp);
    %subtract from all connections to obtain br-en connections
    im_enconn = logical(im_skelconn - im_brconn);
    
%     %--debug--
%     subplot(2,2,4);
%     imshow(im_enconn, [], 'InitialMagnification', 'fit');
    
    %look at length of each br-en connection and count longer connections
    %(relative to the total skeleton length)
    regions_enconn = regionprops(im_enconn);
    peakNum = 0;
    skelLength = cell2mat(struct2cell(regionprops(im_sm,'Area')));
    peakThreshold = 0.07;
    for k = 1:length(regions_enconn)
        area_en = regions_enconn(k).Area;
        if area_en/skelLength > peakThreshold
            peakNum = peakNum + 1;
        end
    end
    
    %if shape is a line, count as two arms
    if peakNum == 1, peakNum = 2; end
    
%     %--debug--
%     title(['shape input has ' num2str(peakNum) ' arms'])
    
end

end


