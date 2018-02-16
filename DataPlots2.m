

mean_arms = zeros(1,3);
std_arms = zeros(1,3);


for clipNum = 1:3
    
    eval(['load mat/Clip' num2str(clipNum) 'Data.mat']);
    eval(['load mat/av_dsp_' num2str(clipNum) '.mat']);
    
    eval(['cent_hist = dataOutput_' num2str(clipNum) '.centroids;']);
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
    
    
    
    eval(['arm_hist = dataOutput_' num2str(clipNum) '.arms;']);
    arms = arm_hist{ind(1)}(1:dataPoints,2);
    for i = 2:length(ind)
        arms = [arms; arm_hist{ind(i)}(1:dataPoints,2)];
    end
    mean_arms(clipNum) = mean(arms);
    std_arms(clipNum) = std(arms);
    
    
    
end

figure;
hold on;
bar(1,mean_arms(1), 'FaceColor', 'k', 'FaceAlpha', 0.6);
bar(2,mean_arms(2), 'FaceColor', 'b', 'FaceAlpha', 0.6);
bar(3,mean_arms(3), 'FaceColor', 'r', 'FaceAlpha', 0.6);
errorbar(mean_arms,std_arms, 'k.')
hold off;
title('Mean no. of arms');
xticks([1 2 3]);
xticklabels({'Control','Chemotaxis','Electrotaxis'})

SaveAsPngEpsandFig(-1,'plots/04', 10, 7/5, 12);

figure; 
loglog(0:length(av_dsp_1)-1, av_dsp_1, 'k', ...
    0:length(av_dsp_2)-1, av_dsp_2, 'b',...
    11:length(av_dsp_3)-1, av_dsp_3(12:end), 'r',...
    0:12, av_dsp_3(1:13), 'r--'); 
legend('Control', 'Chemotaxis', 'Electrotaxis', 'Location', 'best'); 
xlabel('$$log(t)$$', 'Interpreter', 'latex'); ylabel('$$log(\langle x^2 \rangle)$$', 'Interpreter', 'latex');
title('log-log plot of MSD');

SaveAsPngEpsandFig(-1,'plots/05', 10, 7/5, 12);



% arm_hist = dataOutput.arms;
% 
% arms = arm_hist{ind(1)}(1:dataPoints,2);
% 
% figure;
% plot(1:dataPoints, arms, '--');
% text(dataPoints, arms(end), num2str(ind(1)));
% hold on;
% av_arms = arms;
% 
% for i = 2:length(ind)
%     arms = arm_hist{ind(i)}(1:dataPoints,2);
%     
%     plot(1:length(arms), arms, '--');
%     text(dataPoints, arms(end), num2str(ind(i)));
%     
%     av_arms = ((i-1)*av_arms + arms)/(i);
%     
% end
% 
% plot(1:dataPoints, av_arms, 'k');
% text(dataPoints, av_arms(end), 'avg');
% title('Number of arms');
% hold off;
% 
% mean_arms = mean(av_arms);
% std_arms = std(av_arms);
% 
% bar(1:1,mean_arms)
% hold on;
% errorbar(1:1,mean_arms,std_arms)
% hold off;