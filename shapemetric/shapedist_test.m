%load test image
im = im2bw(imread('test.png'));

%show binary image
figure(1);
set(gcf, 'Position', get(0,'Screensize'));
subplot(2,2,1)
imshow(im, [], 'InitialMagnification', 'fit');
title('Original binary image')

%test function    
im_sd = shapedist(im);

%show output in greyscale
subplot(2,2,2)
imshow(mat2gray(im_sd), [], 'InitialMagnification', 'fit');
title('Greyscale shapedist');

%extract boundary of output and show
subplot(2,2,3)
imshow(im_sd.*int8(bwperim(im)), [], 'InitialMagnification', 'fit');
title('Boundary values');

boundary = bwboundaries(im);
boundary = cell2mat(boundary);
ind = sub2ind(size(im), boundary(:,1), boundary(:,2));

x = linspace(0,1,length(ind));
y = double(im_sd(ind))/double(max(im_sd(ind)));

%rearrange plot data y, start at minimum value
tmp = find(y == min(y));
y = [y(tmp:end); y(1:tmp-1)];

%plot values
subplot(2,2,4);
plot(x, y); 


%find and reject peaks with low prominence
peakThreshold = 0.15;
[~,~,w,p] = findpeaks(y,x,'Annotate','extents');
p = p(p > peakThreshold);
peaksNum = length(p);

title(['shape input has ' num2str(peaksNum) ' arms']);









