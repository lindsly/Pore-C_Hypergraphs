%% Load datasets
addpath(genpath('G:\My Drive\UMICH\Research\PoreC\MATLAB\Public_RNA_seq'))
addpath('G:\My Drive\UMICH\Research\PoreC\MATLAB')

% type 1 is IMR-90 and IR fibroblast, type 2 is GM12878

% load('edgeTable_MapQgeq40.mat')

if (exist('Btable','var')) == 0
    load('HumanTfB.mat')
end

if exist('gm12878_RNA_rep1','var') == 0
    load('GM12878_IMR90_RNA.mat');
end

load('samp_nums_TFactories.mat')

tic
for type = 1:2
    if type == 1 && exist('edgeTable','var') == 0
        load('edgeTable_Fibroblast_Combined_enhancer.mat')
    elseif type == 2 && exist('edgeTable_GM12878','var') == 0
    %     load('edgeTable_GM12878.mat')
        load('edgeTable_GM12878_new.mat','edgeTable_GM12878');
    end

    % load('Fibroblast_Adult_Skin_19916.mat') 
    % load('.mat') %% GM12878 RNA-seq FPKM

    %% Format Gene expression
    % Load gene table data
    geneLoc = readtable('biomart_hg38_protein_coding.txt', ...
        'ReadVariableNames', true, 'PreserveVariableNames', true);
    geneMatrix = table2array(geneLoc(:, [2 3 4]));
    allGenes = table2cell(geneLoc(:, 1));

    if type == 1
        rna_table1 = imr90_RNA_rep1(ismember(imr90_RNA_rep1.gene_id,allGenes),:);
        rna_table2 = imr90_RNA_rep2(ismember(imr90_RNA_rep2.gene_id,allGenes),:);
    elseif type == 2
        rna_table1 = gm12878_RNA_rep1(ismember(gm12878_RNA_rep1.gene_id,allGenes),:);
        rna_table2 = gm12878_RNA_rep2(ismember(gm12878_RNA_rep2.gene_id,allGenes),:);
    end

    totalGenes = rna_table1.gene_id;
    geneTPM = mean([rna_table1.TPM,rna_table2.TPM],2);

    %%
    % totalGenes = Fibroblast_Adult_Skin_19916.Properties.RowNames;
    totalGenesTF = Btable.Properties.RowNames;
    totalTF = Btable.Properties.VariableNames;
    % geneTPM = Fibroblast_Adult_Skin_19916.FPKM;
    % edgeTable_w_genes_TF = edgeTable;

    tic
    for p = 1:1000
        p
        if type == 1
            edgeTable_w_genes_TF = edgeTable;
        elseif type == 2
            edgeTable_w_genes_TF = edgeTable_GM12878;
        end


        % % Filter by atac-seq (all nonzero)
        % atacIdx = cellfun(@(x) all(x>0), edgeTable_w_genes_TF.('atac-seq'));
        % edgeTable_w_genes_TF = edgeTable_w_genes_TF(atacIdx, :);
        % disp(['HOCs with accessible chromatin: ', num2str(size(edgeTable_w_genes_TF,1))])
        % 
        % % Filter by chip-seq (at least one nonzero)
        % chipIdx = cellfun(@(x) sum(x)>0, edgeTable_w_genes_TF.('chip-seq'));
        % edgeTable_w_genes_TF = edgeTable_w_genes_TF(chipIdx, :);
        % disp(['HOCs with >= 1 locus with Pol2 binding: ', num2str(size(edgeTable_w_genes_TF,1))])


        %% Get correct number of 3, 4, 5, and 6+ order contacts for sampling
        if type == 1
            samp_num = samp_num_fib;
        elseif type == 2
            samp_num = samp_num_gm;
        end

        samp_idx = [];
        for i_order = 3:length(samp_num)
            idx_temp = find(cellfun('length', edgeTable_w_genes_TF.chrNum) == i_order);
            samp_idx = [samp_idx; datasample(idx_temp, samp_num(i_order),'Replace',false)];
        end
        edgeTable_w_genes_TF = edgeTable_w_genes_TF(samp_idx,:);

        % Check gene expression and common TF
        edgeGenes = edgeTable_w_genes_TF.genes;
        geneExpress = cell(length(edgeGenes), 1);
        comTF = cell(length(edgeGenes), 1);
        enhancer_test = edgeTable_w_genes_TF.enhancer;
        parfor i = 1:length(edgeGenes)
            genes = cellfun(@num2str,edgeGenes{i},'UniformOutput',false);
            geneIdx = ismember(totalGenes, genes);
            geneExpress{i} = geneTPM(geneIdx);
            if length(genes)>=2
                geneTFIdx = ismember(totalGenesTF, genes);
                tfBind = table2array(Btable(geneTFIdx, :))>=3;
                tfIdx = sum(tfBind, 1)>=length(genes);
                comTFTotal = totalTF(tfIdx)';
                tfExpress = geneTPM(ismember(totalGenes, comTFTotal));
                tfExpressIdx = tfExpress>=1;
                comTF{i} = comTFTotal(tfExpressIdx);
            end
        end
        toc

        edgeTable_w_genes_TF.comTF = comTF;

        edgeTable_w_genes_TF.geneExpress = geneExpress;
        num1GenesGT1 = cellfun(@(x) sum(x>=1)>=1, edgeTable_w_genes_TF.geneExpress);
        edgeTable_w_1genes_TF = edgeTable_w_genes_TF(num1GenesGT1,:);

        num2GenesGT1 = cellfun(@(x) sum(x>=1)>=2, edgeTable_w_genes_TF.geneExpress);

        num_enhancer = cellfun(@(x) sum(x>=1)>=1, edgeTable_w_genes_TF.enhancer);

        edgeTable_w_2genes_TF = edgeTable_w_genes_TF(num2GenesGT1,:);
        
        edgeTable_w_2genes_enhancer_TF = edgeTable_w_genes_TF(num2GenesGT1&num_enhancer,:);


        numGenesComTF = cellfun(@(x) size(x,1)>0,edgeTable_w_2genes_enhancer_TF.comTF);

        edgeTable_w_2genes_TF_filter = edgeTable_w_2genes_enhancer_TF(numGenesComTF,:);

        %% Extracing Master Regulators
        tft=Btable(Btable.Properties.VariableNames,Btable.Properties.VariableNames);
        tft{:,:}=(tft{:,:}>=3);
        mr=tft.Properties.RowNames(logical(diag(tft{:,:})));

        tfs_in_HOC = edgeTable_w_2genes_TF_filter.comTF;

        temp = cellfun(@(x) ismember(x,mr),tfs_in_HOC,'UniformOutput',false);
        temp2 = logical(cellfun(@sum, temp)>0);

        edgeTable_w_2genes_TF_filter_MRs = edgeTable_w_2genes_TF_filter(temp2,:);

        %% Extracting Gene Table
        if type == 1
            HOCs_order = cellfun('length', edgeTable.chrNum);
        elseif type == 2
            HOCs_order = cellfun('length',edgeTable_GM12878.chrNum);
        end
        Tfactories_order = cellfun('length', edgeTable_w_genes_TF.chrNum);
        gt1genes_order = cellfun('length', edgeTable_w_1genes_TF.chrNum);
        gt2genes_order = cellfun('length', edgeTable_w_2genes_TF.chrNum);
        commonTFs_order = cellfun('length', edgeTable_w_2genes_TF_filter.chrNum);
        MRs_order = cellfun('length', edgeTable_w_2genes_TF_filter_MRs.chrNum);
        enhancer_order = cellfun('length',edgeTable_w_2genes_enhancer_TF.chrNum);

        HOC_table{p} = [cntr(HOCs_order==3) cntr(Tfactories_order==3) cntr(gt2genes_order==3) cntr(enhancer_order==3) cntr(commonTFs_order==3) cntr(MRs_order==3);...
                     cntr(HOCs_order==4) cntr(Tfactories_order==4) cntr(gt2genes_order==4) cntr(enhancer_order==4) cntr(commonTFs_order==4) cntr(MRs_order==4);...
                     cntr(HOCs_order==5) cntr(Tfactories_order==5) cntr(gt2genes_order==5) cntr(enhancer_order==5) cntr(commonTFs_order==5) cntr(MRs_order==5);...
                     cntr(HOCs_order>=6) cntr(Tfactories_order>=6) cntr(gt2genes_order>=6) cntr(enhancer_order>=6) cntr(commonTFs_order>=6) cntr(MRs_order>=6)];

    end

    if type == 1
        HOC_table_Fib = HOC_table;
    elseif type == 2
        HOC_table_GM = HOC_table;
    end
    toc
end


save('HOC_table_GM_1000trials.mat','HOC_table_GM');
save('HOC_table_Fib_1000trials.mat','HOC_table_Fib');

