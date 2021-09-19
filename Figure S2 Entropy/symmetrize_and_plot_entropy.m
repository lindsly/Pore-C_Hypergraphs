%% Load entropy data and symmetrize
load('fibroblast_1Mb_3_way_entropy.mat')
sym_fib = double(full(symtensor(tensor(entropyChr./log(nodeSize)))));
fib_intra_entropy = [];
for i = 1:24
    fib_intra_entropy = [fib_intra_entropy; entropyChr(i,i,i)];
end

load('GM12878_1Mb_3_way_entropy.mat')
sym_GM = double(full(symtensor(tensor(entropyChr./log(nodeSize)))));
gm_intra_entropy = [];
for i = 1:23
    gm_intra_entropy = [gm_intra_entropy; entropyChr(i,i,i)];
end


%% Plotting
red_map = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
clim_vals = [0 max(max(max(max(sym_fib))),max(max(max(sym_GM))))];

figure
subplot(1,2,1)
vol3d('CData',sym_fib), view(3)
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Box','on','BoxStyle', ...
    'full','XLim',[0 24],'YLim',[0 24],'ZLim',[0 24],'CLim',clim_vals)
axis square
title('Fibroblast Entropy')
subplot(1,2,2)
vol3d('CData',sym_GM), view(3)
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Box','on','BoxStyle', ...
    'full','XLim',[0 23],'YLim',[0 23],'ZLim',[0 23],'CLim',clim_vals)
axis square
colormap(red_map)
title('B Lymphocyte Entropy')

%% Intra-chr entropy
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
             '17' '18' '19' '20' '21' '22' 'X' 'Y'}';

%% Upside down     
figure('Position', [540 288 771 420])
hold on
b1 = bar(fib_intra_entropy,'BarWidth', .6,'FaceColor','flat');
b2 = bar(-[gm_intra_entropy; 0],'BarWidth', .6,'FaceColor','flat');

b1.CData = ten_colors(9,:);
b2.CData = ten_colors(8,:);

ylim_temp = get(gca,'YLim');
ylim_temp(1) = ylim_temp(1)-.1;
ylim_temp(2) = ylim_temp(2)+.1;

set(gca,'YLim',ylim_temp,'TickLength',[0 0],'XTick',1:24,...
    'XTickLabel',chr_nums,'LineWidth',.1,'FontSize',14)
xlabel('Chromosome')
box on

%% Side by side  
close all
figure('Position', [195 288 1509 420])
hold on
b = bar([fib_intra_entropy,[gm_intra_entropy; 0]],'BarWidth', .9,'FaceColor','flat');

b(1).CData = ten_colors(8,:);
b(2).CData = ten_colors(9,:);

ylim_temp = get(gca,'YLim');
ylim_temp(1) = 2;
ylim_temp(2) = ylim_temp(2)+.1;

set(gca,'YLim',ylim_temp,'TickLength',[0 0],'XTick',1:24,...
    'XTickLabel',chr_nums,'LineWidth',.1,'FontSize',14)
xlabel('Chromosome')
box on