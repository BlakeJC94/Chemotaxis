%load test image
im = im2bw(imread('test1.png'));

%show binary image
figure(1);
% set(gcf, 'Position', get(0,'Screensize'));
% subplot(2,2,1)
imshow(im, [], 'InitialMagnification', 'fit');
title('Original binary image')

%test function
numPeaks = shapemetric(im);

title(['shape input has ' num2str(numPeaks) ' arms']);
