load('degree_chosen_25mb_both_cell_types.mat')

% chrLength = ceil([249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
%     134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58]./25);
% chrLength_afterzeros = cumsum(chrLength) - [1 1 1 1 2 2 2 2 2 2 2 3 3 3 4 4 4 5 5 5 5 6 6 8];

% Remove indices that are zero in both
% idx = (degree_chosen_fib_25mb+degree_chosen_gm_25mb)==0;
% degree_chosen_fib_25mb(idx) = [];
% degree_chosen_gm_25mb(idx) = [];

% Remove indices that are zero in each separately
degree_chosen_fib_25mb(degree_chosen_fib_25mb ==0) = [];
degree_chosen_gm_25mb(degree_chosen_gm_25mb ==0) = [];

% top_val = max(max(max(degree_chosen_fib_25mb)),max(max(degree_chosen_gm_25mb)));

degree_chosen_fib_25mb=degree_chosen_fib_25mb/max(degree_chosen_fib_25mb);
degree_chosen_gm_25mb=degree_chosen_gm_25mb/max(degree_chosen_gm_25mb);

figure('Position', [395 504 1197 420])
bar([degree_chosen_fib_25mb],'BarWidth',.6,'FaceColor','k')
set(gca,'TickLength',[0 0],'YLim',[0 1.1])

figure('Position', [395 504 1197 420])
bar([degree_chosen_gm_25mb],'BarWidth',.6,'FaceColor','k') %[225/255,110/255,40/255]
set(gca,'TickLength',[0 0],'YLim',[0 1.1])