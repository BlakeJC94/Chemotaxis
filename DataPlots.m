% if ~exist('dataOutput','var')
%     dataOutput = AnalyseDictyShape();
% end

clipNum = 1;

eval(['load mat/Clip' num2str(clipNum) 'Data.mat']);
eval(['dataOutput = dataOutput_' num2str(clipNum)]);




cent_hist = dataOutput.centroids;
FRAME_RANGE = dataOutput.Frames;
imsize = dataOutput.imsize;

close all;
 

%% Filter out cells that have useful data
% Reject cells that have only been recorded for less than half of the video
% (minPoints) and reject cells that have moved less than 10 pixels from 
% initial position by final frame


thrDisp = 30;

% dataPoints = floor(numel(FRAME_RANGE(1):FRAME_RANGE(2))/2);
dataPoints = 40;
ind = zeros(1,size(cent_hist,2));
for i = 1:size(cent_hist,2)
    
    totalDisplacement = norm(cent_hist{i}(end,2:3) - cent_hist{i}(1,2:3));
    if (size(cent_hist{i},1) >= dataPoints) && (totalDisplacement > thrDisp)
        ind(i) = 1;
    end
    
end

ind = find(ind);

dataPoints = size(cent_hist{ind(1)}, 1);
for i = 2:length(ind)
    dataPoints = min(dataPoints, size(cent_hist{ind(i)}, 1));
end



%% Plot displacement trajectories
figure;

x = cent_hist{ind(1)}(:,2) - cent_hist{ind(1)}(1,2);
y = - (cent_hist{ind(1)}(:,3) - cent_hist{ind(1)}(1,3));

plot(x,y);
% text(x(end), y(end), num2str(ind(1)));

xlim = max(abs(x));
ylim = max(abs(y));

hold on;

for i = 2:length(ind)
    
    x = cent_hist{ind(i)}(:,2) - cent_hist{ind(i)}(1,2);
    y = - (cent_hist{ind(i)}(:,3) - cent_hist{ind(i)}(1,3));
    
    xlim = max(xlim, max(abs(x)));
    ylim = max(ylim, max(abs(y)));
    
    plot(x,y);
%     text(x(end), y(end), num2str(ind(i)));
    
end
hold off;
title('Cell trajectories');
axis(1.1*[-xlim xlim -ylim ylim]);
xlabel('$$x$$', 'Interpreter', 'latex'); ylabel('$$y$$', 'Interpreter', 'latex'); 

SaveAsPngEpsandFig(-1,'plots/01', 7, 7/5, 12);


%% Plot Directionality
figure;

x = linspace(-1,1,100);
y1 = +sqrt(1 - x.^2);
y2 = -sqrt(1 - x.^2);

plot(x, y1, 'k', x, y2, 'k');
axis square;
hold on;

vec = [cent_hist{ind(1)}(end,2) - cent_hist{ind(1)}(1,2),...
     - (cent_hist{ind(1)}(end,3) - cent_hist{ind(1)}(1,3))];
vec = vec/norm(vec);
    
plot([0, vec(1)], [0, vec(2)], '--');
% text(vec(1), vec(2), num2str(ind(1)));

av_vec = vec;

for i = 2:length(ind)
    
    vec = [cent_hist{ind(i)}(end,2) - cent_hist{ind(i)}(1,2),...
         - (cent_hist{ind(i)}(end,3) - cent_hist{ind(i)}(1,3))];
    vec = vec/norm(vec);
    
    plot([0, vec(1)], [0, vec(2)], '--');
%     text(vec(1), vec(2), num2str(ind(i)));
    
    av_vec = ((i-1)*av_vec + vec)/(i);
    
end

quiver(0,0,av_vec(1), av_vec(2), 'k', 'LineWidth', 7)
% plot([0 av_vec(1)], [0 av_vec(2)], 'k' );


hold off;
title('Cell directionality');

xlabel('$$x$$', 'Interpreter', 'latex'); ylabel('$$y$$', 'Interpreter', 'latex'); 

SaveAsPngEpsandFig(-1,'plots/02', 5, 1, 12);



%% Plot RMS displacement
figure; 

dsp = [cent_hist{ind(1)}(:,2) - cent_hist{ind(1)}(1,2), ...
    (cent_hist{ind(1)}(:,3) - cent_hist{ind(1)}(1,3))];
dsp = sum(dsp(1:dataPoints,:).^2,2);

plot(1:dataPoints, dsp, '--');
% text(dataPoints, dsp(end), num2str(ind(1)));
hold on; 

av_dsp = dsp;

for i = 2:length(ind)
    dsp = [cent_hist{ind(i)}(:,2) - cent_hist{ind(i)}(1,2), ...
        (cent_hist{ind(i)}(:,3) - cent_hist{ind(i)}(1,3))];
    dsp = sum(dsp(1:dataPoints, :).^2,2);
    
    plot(0:length(dsp)-1, dsp, '--');
%     text(dataPoints, dsp(end), num2str(ind(i)));
    
    av_dsp = ((i-1)*av_dsp + dsp)/(i);
    
end


plot(0:dataPoints-1, av_dsp, 'k');
text(dataPoints, av_dsp(end), 'avg');

hold off; 
% title(['Mean squared displacement over ' num2str(dataPoints) ' frames']);
title(['Squared displacement']);
xlabel('$$t$$ (frames)', 'Interpreter', 'latex'); ylabel('pixels$$^2$$', 'Interpreter', 'latex'); 

SaveAsPngEpsandFig(-1,'plots/03', 7, 7/5, 12);

% figure; 
% plot(1:dataPoints, sqrt(av_dsp), 'k');
% title(['Root mean displacement over ' num2str(dataPoints) ' frames']);


%% Plot Mean number of arms??

arm_hist = dataOutput.arms;

arms = arm_hist{ind(1)}(1:dataPoints,2);

figure;
plot(1:dataPoints, arms, '--');
text(dataPoints, arms(end), num2str(ind(1)));
hold on;
av_arms = arms;

for i = 2:length(ind)
    arms = arm_hist{ind(i)}(1:dataPoints,2);
    
    plot(1:length(arms), arms, '--');
    text(dataPoints, arms(end), num2str(ind(i)));
    
    av_arms = ((i-1)*av_arms + arms)/(i);
    
end

plot(1:dataPoints, av_arms, 'k');
text(dataPoints, av_arms(end), 'avg');
title('Number of arms');
hold off;

mean_arms = mean(av_arms);
std_arms = std(av_arms);

bar(1:1,mean_arms)
hold on;
errorbar(1:1,mean_arms,std_arms)
hold off;

