%% Load Data
load('comb1234_chr22_100kb.mat')
red_map = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];

% [hic,bad_locs] = hic_trim(comb1234_chr22_100kb,1,.1);

%% TAD Calculation
% Calculate TADs via Jie et al. 2016
TAD_boundaries = TAD_Laplace(comb1234_chr22_100kb, 0.8, 3, 1, 1);

% Adjust for centromere
TAD_boundaries = TAD_boundaries+105;

% Adjust for region of interest
TAD_boundaries = TAD_boundaries(60:62)-333;

%% Plotting
figure
imagesc(log(comb1234_chr22_100kb(229:241, 229:241))) 
plot_TADs(TAD_boundaries)
axis square
set(gca,'TickLength',[0 0],'XTick',1:13,'XTickLabel',...
        229:241,'YTick',1:13,'YTickLabel',229:241)
colormap(red_map)

% %% temp plotting
% figure
% % imagesc(log(comb1234_chr22_100kb+.5)+1)
% imagesc((comb1234_chr22_100kb))
% colormap(red_map)
% plot_TADs(TAD_boundaries)
% axis square
% 
% %% temp plotting2
% figure
% imagesc(log(hic+.5)+1)
% colormap(red_map)
% plot_TADs(TAD_boundaries)
% axis square
