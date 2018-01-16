%load test image
im = im2bw(imread('test7.png'));
im = imfill(im, 'holes');

%show binary image
figure(1);
set(gcf, 'Position', get(0,'Screensize'));
subplot(1,2,1)
imshow(im, [], 'InitialMagnification', 'fit');
title('Original binary image')

%test function
[peakNum, im_enconn] = skelmetric(im);

subplot(1,2,2)
imshow(im_enconn, [], 'InitialMagnification', 'fit');


text(0,0,['end-branch connections, ' num2str(peakNum) ' arms'],'Color','g');









