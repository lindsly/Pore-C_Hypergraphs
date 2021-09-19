%% Load distance data
load('chrDistance.mat')
load('nullDistri.mat')

figure
subplot(1,2,1)
% Stem plot
hist(nullDistri)
title('Null Distribution')
subplot(1,2,2)
bar(chrDistance)
title('Chromosome Distances')


%% Chr Dist plots
ten_colors = [       0    0.4470    0.7410; % 1
                0.8500    0.3250    0.0980; % 2
                0.9290    0.6940    0.1250; % 3
                0.4940    0.1840    0.5560; % 4
                0.4660    0.6740    0.1880; % 5
                0.3010    0.7450    0.9330; % 6
                0.6350    0.0780    0.1840; % 7
                0.9000    0.4000    0.4000; % 8
                0.4000    0.4000    0.8000; % 9
                0.5000    0.500     0.5000; % 10
                     0         0         0;]; % 11
 chr_nums = {'1' '2' '3' '4' '5' '6' '7' '8' '9'...
             '10' '11' '12' '13' '14' '15' '16' ...
             '17' '18' '19' '20' '21' '22' 'X'}';

figure('Position', [564 269 771 420])
hold on
b = bar(chrDistance,'BarWidth', .6,'FaceColor','flat');

b.CData = ten_colors(5,:);

ylim_temp = get(gca,'YLim');
% ylim_temp(1) = 2;
% ylim_temp(2) = ylim_temp(2)+.1;

set(gca,'YLim',ylim_temp,'TickLength',[0 0],'XTick',1:23,...
    'XTickLabel',chr_nums,'LineWidth',.1,'FontSize',14)
xlabel('Chromosome')
box on

% p_vals = sum(chrDistance(i)<=nullDistri)/length(nullDistri);