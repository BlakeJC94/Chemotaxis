function [numPeaks, im_sd] = shapemetric(im)
% Performs shapedist transform on binary image of single region and
% identifies number of arms a simply connected irregular polygon has.
%
% Function will accept binary image of single region and return the number
% of arms, based on prominence of peaks in values of shapedist on
% boundaries
%
% INPUTS:
%   im - binary image of a region
%
% OUTPUTS:
%   numPeaks - number of peaks in plot of boundary vs shapedist
%
% Blake Cook 2018-01-04
% ------------------------

%get shapedist transform
im_sd = shapedist(im);

%get indicies of image boundary values
boundary = bwboundaries(im);
boundary = cell2mat(boundary);
ind = sub2ind(size(im), boundary(:,1), boundary(:,2));

%get plot values, extract boundaries from im_sd
x = linspace(0,1,length(ind));
y = double(im_sd(ind))/double(max(im_sd(ind)));

%rearrange plot data y, start at minimum value
tmp = find(y == min(y));
y = [y(tmp:end); y(1:tmp-1)];

% %-- DEBUG: plot values
% subplot(2,2,4);
% plot(x, y); 

%find and reject peaks with low prominence
peakThreshold = 0.15;
[~,~,~,p] = findpeaks(y,x);
p = p(p > peakThreshold);
numPeaks = length(p);


end


