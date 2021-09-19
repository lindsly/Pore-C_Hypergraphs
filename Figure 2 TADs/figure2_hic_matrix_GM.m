%% Load Data
load('GM12878_chr22_100kb.mat')
red_map = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];

%% TAD Calculation
% Calculate TADs via Jie et al. 2016
TAD_boundaries = TAD_Laplace(GM12878_chr22_100kb, 0.8, 3, 1, 1);

% Adjust for centromere
TAD_boundaries = TAD_boundaries+105;

% Adjust for region of interest
TAD_boundaries = TAD_boundaries(5:7)-117;

%% Plotting
figure
imagesc(log(GM12878_chr22_100kb))%(13:22, 13:22)
plot_TADs(TAD_boundaries)
axis square
set(gca,'TickLength',[0 0],'XTick',1:10,'XTickLabel',...
        118:127,'YTick',1:10,'YTickLabel',118:127)
colormap(red_map)