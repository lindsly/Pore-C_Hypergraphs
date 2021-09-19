%% Load data and parameters

% Load interChrContacts
load('gmInterChrContacts25mb.mat')

chrLength = ceil([249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
            134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58]'/25);
        
% sum(chrLength)
        
inter_opts = {'','chrs_2','chrs_3','chrs_4','chrs_5'};
inter_names = {'','Chromosomes_2','Chromosomes_3','Chromosomes_4','Chromosomes_5'};

%% Find weights of all inter-chr interactions per chromosome
interSize = cellfun('length', interChrContacts.interChrSizeUnique);
for iInter = 2:5
    interChr.(inter_opts{iInter}) = interChrContacts(interSize == iInter,:);
    chr_weights.(inter_opts{iInter}) = zeros(24,1);
    for iChr = 1:24
        chr_locs = cellfun(@(x) ismember(x,iChr), interChr.(inter_opts{iInter}).interChrSizeUnique, 'UniformOutput', false);
        idx = logical(cellfun(@(x) sum(x), chr_locs));
        chr_weights.(inter_opts{iInter})(iChr,1) = sum(interChr.(inter_opts{iInter}).weights(idx));
    end
end

%% Plot Chromosome Order
figure('Position', [11 433 1850 420])
subplot(1,4,1)
    plot(chr_weights.chrs_2,'.','MarkerSize',10)
    title('2 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,2)
    plot(chr_weights.chrs_3,'.','MarkerSize',10)
    title('3 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,3)
    plot(chr_weights.chrs_4,'.','MarkerSize',10)
    title('4 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,4)
    plot(chr_weights.chrs_4,'.','MarkerSize',10)
    title('4 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,4)
    plot(chr_weights.chrs_5,'.','MarkerSize',10)
    title('5 Chromosomes')
    xlabel('Chromosome')
suptitle('Chromosome Order')

%% Plot Chr Length Order
% [~,idx] = sort(chrLength,'descend');
% 
% figure('Position', [11 433 1850 420])
% subplot(1,3,1)
%     plot(chr_weights.chrs_2(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('2 Chromosomes')
%     xlabel('Chromosome')
% subplot(1,3,2)
%     plot(chr_weights.chrs_3(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('3 Chromosomes')
%     xlabel('Chromosome')
% subplot(1,3,3)
%     plot(chr_weights.chrs_4(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('4 Chromosomes')
%     xlabel('Chromosome')
%     
% suptitle('Chromosome Length Order')

%% Reduction of inter-chr interactions by prevalance of chromosomes
scaling_factors = [0 1E5 5E9 1E14]
interChr_norm = interChr;
for iInter = 2:5
    interChr_norm.(inter_opts{iInter}).weights_norm = interChr_norm.(inter_opts{iInter}).weights;
    for iChr = 1:24
        chr_locs = cellfun(@(x) ismember(x,iChr), interChr.(inter_opts{iInter}).interChrSizeUnique, 'UniformOutput', false);
        idx = logical(cellfun(@(x) sum(x), chr_locs));
        % Divide by total contacts involving each chromosome
%         interChr_norm.(inter_opts{iInter}).weights_norm(idx) = ...
%             interChr_norm.(inter_opts{iInter}).weights_norm(idx)...
%             /chr_weights.(inter_opts{iInter})(iChr);
        % Divide by chromosome length
        interChr_norm.(inter_opts{iInter}).weights_norm(idx) = ...
            interChr_norm.(inter_opts{iInter}).weights_norm(idx)...
            /chrLength(iChr);
    end
%     interChr_norm.(inter_opts{iInter}).weights_norm =...
%         interChr_norm.(inter_opts{iInter}).weights_norm*scaling_factors(iInter);
end

%% Find weights of all inter-chr interactions per chromosome after normalization
for iInter = 2:5
    for iChr = 1:24
        chr_locs = cellfun(@(x) ismember(x,iChr), interChr_norm.(inter_opts{iInter}).interChrSizeUnique, 'UniformOutput', false);
        idx = logical(cellfun(@(x) sum(x), chr_locs));
        chr_weights_norm.(inter_opts{iInter})(iChr,1) = sum(interChr_norm.(inter_opts{iInter}).weights_norm(idx));
    end
end

%% Plot Chromosome Order after Normalization
figure('Position', [11 433 1850 420])
subplot(1,4,1)
    plot(chr_weights_norm.chrs_2,'.','MarkerSize',10)
    title('2 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,2)
    plot(chr_weights_norm.chrs_3,'.','MarkerSize',10)
    title('3 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,3)
    plot(chr_weights_norm.chrs_4,'.','MarkerSize',10)
    title('4 Chromosomes')
    xlabel('Chromosome')
subplot(1,4,4)
    plot(chr_weights_norm.chrs_5,'.','MarkerSize',10)
    title('5 Chromosomes')
    xlabel('Chromosome')
suptitle('Chromosome Order')

%% Plot Chr Length Order after Normalization
% [~,idx] = sort(chrLength,'descend');
% 
% figure('Position', [11 433 1850 420])
% subplot(1,4,1)
%     plot(chr_weights_norm.chrs_2(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('2 Chromosomes')
%     xlabel('Chromosome')
% subplot(1,4,2)
%     plot(chr_weights_norm.chrs_3(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('3 Chromosomes')
%     xlabel('Chromosome')
% subplot(1,4,3)
%     plot(chr_weights_norm.chrs_4(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('4 Chromosomes')
%     xlabel('Chromosome')
% subplot(1,4,4)
%     plot(chr_weights_norm.chrs_5(idx),'.','MarkerSize',10)
%     xticks(1:24)
%     xticklabels(idx)
%     title('5 Chromosomes')
%     xlabel('Chromosome')
%     
% % suptitle('Chromosome Length Order - Divide by Total Weight')
% suptitle('Chromosome Length Order - Divide by Chr Length')

%%
all_chrs_max = cell(72,1);
all_chrs_max_table = table;
all_chrs_norm_weight = zeros(72,1);
all_chrs_max_table.Chromosome = [1:23]';
counter = 1;
for iInter = 2:5
    chr_specific_inter.(inter_opts{iInter}) = cell(23,1);
    all_chrs_max_temp = cell(23,1);
    for iChr = 1:23
        chr_locs = cellfun(@(x) ismember(x,iChr), interChr_norm.(inter_opts{iInter}).interChrSizeUnique, 'UniformOutput', false);
        idx = logical(cellfun(@(x) sum(x), chr_locs));
        chr_specific_inter.(inter_opts{iInter}){iChr} = interChr_norm.(inter_opts{iInter})(idx,:);
        
        [~,max_idx] = max(chr_specific_inter.(inter_opts{iInter}){iChr}.weights_norm);
        all_chrs_max(counter) = chr_specific_inter.(inter_opts{iInter}){iChr}.interChrSizeUnique(max_idx);
        all_chrs_norm_weight(counter) = chr_specific_inter.(inter_opts{iInter}){iChr}.weights_norm(max_idx);
%         all_chrs_weight(counter) = chr_specific_inter.(inter_opts{iInter}){iChr}.we
        counter = counter + 1;
        all_chrs_max_temp(iChr) = chr_specific_inter.(inter_opts{iInter}){iChr}.interChrSizeUnique(max_idx);
    end
    all_chrs_max_table.(inter_names{iInter}) = all_chrs_max_temp;
end
temp = cellfun(@(x) sort(char(x)),all_chrs_max,'UniformOutput',false);
[temp, IA] = unique(temp,'rows')

all_chrs_max = all_chrs_max(IA);

all_chrs_max_table_gm = all_chrs_max_table;

% temp = cellstr(all_chrs_max_table.Chromosomes_2)
writetable(all_chrs_max_table, 'Inter_Chromosomal_Interactions_GM.csv');
%% All chromosomes table
