%load test image
im = im2bw(imread('test4.png'));

%test function    
im_sd = shapedist(im);

% disp(['---- shapedist output ----']);
% out = [0:14; (1:14)' im_sd];
% disp(out);

subplot(1,2,1)
imshow(im, [], 'InitialMagnification', 'fit');
title('Original binary image')

subplot(1,2,2)
imshow(mat2gray(im_sd), [], 'InitialMagnification', 'fit');
title('Greyscale shapedist');

figure;
imshow(im_sd.*int8(bwperim(im)), [], 'InitialMagnification', 'fit');

boundary = bwboundaries(im);
boundary = cell2mat(boundary);
ind = sub2ind(size(im), boundary(:,1), boundary(:,2));
im_sd(ind);

plot(1:length(ind), im_sd(ind));