%% Load data
load('final_HOCs_all_filters_Fibroblast.mat')
Fib_TFactory_Common_TFs = edgeTable_w_2genes_TF_filter;
load('final_HOCs_all_filters_GM12878.mat')
GM_TFactory_Common_TFs = edgeTable_w_2genes_TF_filter;
clear edgeTable_w_2genes_TF_filter

chrLength = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, ...
    159345973, 145138636, 138394717, 133797422, 135086622, ...
    133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, ...
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415];

addpath('../Pathways')

if (exist('kegg','var')) == 0
	load('../Pathways/kegg.mat')
    load('../Pathways/kegg_info.mat')
end

pathways_OI = kegg_info.KEGG_Pathway;%{'CellCycle','Circadian','Glycolysis','Krebs','BCR','OxPhos','PentPhos','P53'}; %,'Cancer'
pathways_ID = kegg_info.ID;%{'k04110', 'k04710', 'k00010', 'k00020','k04662','k00190','k00030','k04115'}; %,'k05200'
pathways_OI = regexprep(pathways_OI,{' ','/','-','!',',','(',')'},'_');
pathways_OI = regexprep(pathways_OI,{'\_{2,}'},'_');

Fib_genes = Fib_TFactory_Common_TFs.genes;
GM_genes = GM_TFactory_Common_TFs.genes;

% Load gene table data
geneLoc = readtable('biomart_hg38_protein_coding.txt', ...
    'ReadVariableNames', true, 'PreserveVariableNames', true);
geneLoc.Properties.VariableNames = {'genenames','chr','start','stop'};
for i = 1:size(geneLoc)
    geneLoc.start(i) = geneLoc.start(i) + sum(chrLength(1:geneLoc.chr(i)-1));
    geneLoc.stop(i) = geneLoc.stop(i) + sum(chrLength(1:geneLoc.chr(i)-1));
end

tic
for type = 1:2
    %% Associate Genes from biomart to read IDs
    if type == 1
        genes_info = cell(size(Fib_genes,1),5);
    elseif type == 2
        genes_info = cell(size(GM_genes,1),5);
    end
    
    for i = 1:size(genes_info,1)
        i
        if type == 1
            [C,~,idx] = intersect(Fib_genes{i},geneLoc.genenames,'stable');
        elseif type == 2
            [C,~,idx] = intersect(GM_genes{i},geneLoc.genenames,'stable');
        end
        genes_info{i,1} = C;
        genes_info{i,2} = geneLoc.chr(idx);
        genes_info{i,3} = geneLoc.start(idx);
        genes_info{i,4} = geneLoc.stop(idx);
        genes_info{i,5} = (genes_info{i,4}+genes_info{i,3})/2;

    %     genes_info{i,6} = squareform(pdist(genes_info{i,5}));
    end

    genes_info = cell2table(genes_info,'VariableNames',{'genenames','chr','start','stop','middle'});
    if type == 1
        genes_info.readID = Fib_TFactory_Common_TFs.readID;
    elseif type == 2
        genes_info.readID = GM_TFactory_Common_TFs.readID;
    end
    genes_info = movevars(genes_info,'readID','Before','genenames');

    %% Associate pathways with genes in HOCs
    pathways_temp = cell(size(genes_info,1),1);
    far_apart = zeros(size(genes_info,1),1);
    parfor i = 1:size(genes_info,1)
        i
        [~, LOCB] = ismember(genes_info.genenames{i},kegg.Properties.RowNames);
        LOCB(LOCB == 0) = [];
        if max(sum(kegg{LOCB,:})) > 1
            idx = sum(kegg{LOCB,:}) > 1;
            pathways_temp{i} = kegg_info.KEGG_Pathway(idx);
            far_apart(i) = all(pdist(genes_info.middle{i}) > 1E5);
        end
    end

    genes_info.pathways = pathways_temp;
    genes_info.far_apart = far_apart;
    idx = cellfun('length', pathways_temp) > 0 & far_apart > 0;

    genes_info_OI = genes_info(idx,:);

    %%
    if type == 1
        Fib_TFactory_GOI = Fib_TFactory_Common_TFs(ismember(Fib_TFactory_Common_TFs.readID,genes_info_OI.readID),:);
        Fib_genes_info_OI = genes_info_OI;
    elseif type == 2
        GM_TFactory_GOI = GM_TFactory_Common_TFs(ismember(GM_TFactory_Common_TFs.readID,genes_info_OI.readID),:);
        GM_genes_info_OI = genes_info_OI;
    end
end
toc

%% Current examples
Fib_examples = [6 7 9 10 13 21];
GM_examples = [11 16 29 34 38 44];

Fib_examples = Fib_TFactory_GOI(Fib_examples,:);
GM_examples = GM_TFactory_GOI(GM_examples,:);


%% Ideogram figures for first example

hs_cytobands = cytobandread('hs_cytoBand.txt');
% hfchrom = chromosomeplot_sl(hs_cytobands);
hfchrom = chromosomeplot(hs_cytobands);

set(gcf, 'color', 'none'); set(gca, 'color', 'none');



%% Intersection
% parfor i = 1:size(Fib_genes,1)
%     i
%     gene_int = cellfun(@(x) intersect(Fib_genes{i},x),GM_genes,'UniformOutput',false);
%     idx_geq2(:,i) = cellfun("length", gene_int)>=2;
%     idx_geq3(:,i) = cellfun("length", gene_int)>=3;
%     idx_geq4(:,i) = cellfun("length", gene_int)>=4;
%     idx_geq5(:,i) = cellfun("length", gene_int)>=5;
% end
% 
% [row, col] = find(idx_geq2);
% [row, col] = find(idx_geq3); % No valid results
% [row, col] = find(idx_geq4); % No valid results
