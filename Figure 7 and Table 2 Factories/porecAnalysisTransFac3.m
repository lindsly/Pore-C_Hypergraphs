%% Load datasets
addpath(genpath('G:\My Drive\UMICH\Research\PoreC\MATLAB\Public_RNA_seq'))
addpath('G:\My Drive\UMICH\Research\PoreC\MATLAB')

% type 1 is IMR-90 and IR fibroblast, type 2 is GM12878
type = 1;

% load('edgeTable_MapQgeq40.mat')

if (exist('Btable','var')) == 0
    load('HumanTfB.mat')
end

if exist('gm12878_RNA_rep1','var') == 0
    load('GM12878_IMR90_RNA.mat');
end


% if type == 1 && exist('Fib1','var') == 0
%     load('ENCFF081NBY.mat') % Fibroblast RNA-seq ENCFF956RID is replicate 1 (Fib1)
%     load('ENCFF956RID.mat') % Fibroblast RNA-seq ENCFF081NBY is replicate 2 (Fib2)
% elseif type == 2 && exist('ENCFF306TLL','var') == 0
%     load('ENCFF306TLL.mat') %% GM12878 RNA-seq ENCFF306TLL is replicate 1
%     load('ENCFF418FIT.mat') %% GM12878 RNA-seq ENCFF418FIT is replicate 2
% end

tic
if type == 1 && exist('edgeTable','var') == 0
    load('edgeTable_Fibroblast_Combined.mat')
%     load('edgeTable_Fibroblast_Combined_enhancer.mat')
elseif type == 2 && exist('edgeTable_GM12878','var') == 0
%     load('edgeTable_GM12878.mat')
    load('edgeTable_GM12878_new.mat','edgeTable_GM12878');
end
toc

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
if type == 1
    edgeTable_w_genes_TF = edgeTable;
elseif type == 2
    edgeTable_w_genes_TF = edgeTable_GM12878;
end


% Filter by atac-seq (all nonzero)
atacIdx = cellfun(@(x) all(x>0), edgeTable_w_genes_TF.('atac-seq'));
edgeTable_w_genes_TF = edgeTable_w_genes_TF(atacIdx, :);
disp(['HOCs with accessible chromatin: ', num2str(size(edgeTable_w_genes_TF,1))])

% Filter by chip-seq (at least one nonzero)
chipIdx = cellfun(@(x) sum(x)>0, edgeTable_w_genes_TF.('chip-seq'));
edgeTable_w_genes_TF = edgeTable_w_genes_TF(chipIdx, :);
disp(['HOCs with >= 1 locus with Pol2 binding: ', num2str(size(edgeTable_w_genes_TF,1))])

% Filter by enhancer (at least one nonzero)
% enhancerIdx = cellfun(@(x) sum(x)>0, edgeTable_w_genes_TF.('enhancer'));
% edgeTable_w_genes_TF = edgeTable_w_genes_TF(enhancerIdx, :);
% disp(['HOCs with >= 1 enhancer: ', num2str(size(edgeTable_w_genes_TF,1))])


% Check gene expression and common TF
edgeGenes = edgeTable_w_genes_TF.genes;
geneExpress = cell(length(edgeGenes), 1);
comTF = cell(length(edgeGenes), 1);
parfor i = 1:length(edgeGenes)
%     i
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

edgeTable_w_genes_TF.comTF = comTF;

edgeTable_w_genes_TF.geneExpress = geneExpress;
num1GenesGT1 = cellfun(@(x) sum(x>=1)>=1, edgeTable_w_genes_TF.geneExpress);
edgeTable_w_1genes_TF = edgeTable_w_genes_TF(num1GenesGT1,:);
disp(['HOCs >=1 gene: ', num2str(sum(num1GenesGT1))])

num2GenesGT1 = cellfun(@(x) sum(x>=1)>=2, edgeTable_w_genes_TF.geneExpress);
sum(num2GenesGT1)
disp(['HOCs >=2 gene: ', num2str(sum(num2GenesGT1))])

edgeTable_w_2genes_TF = edgeTable_w_genes_TF(num2GenesGT1,:);
% edgeTable_w_genes_TF = edgeTable_w_genes_TF(numGenesGT1,:);

numGenesComTF = cellfun(@(x) size(x,1)>0,edgeTable_w_2genes_TF.comTF);
disp(['HOCs >=2 gene and common TFs: ', num2str(sum(numGenesComTF))])

edgeTable_w_2genes_TF_filter = edgeTable_w_2genes_TF(numGenesComTF,:);

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


HOC_table_real = [cntr(HOCs_order==3) cntr(Tfactories_order==3) cntr(gt1genes_order==3) cntr(gt2genes_order==3) cntr(commonTFs_order==3) cntr(MRs_order==3);...
             cntr(HOCs_order==4) cntr(Tfactories_order==4) cntr(gt1genes_order==4) cntr(gt2genes_order==4) cntr(commonTFs_order==4) cntr(MRs_order==4);...
             cntr(HOCs_order==5) cntr(Tfactories_order==5) cntr(gt1genes_order==5) cntr(gt2genes_order==5) cntr(commonTFs_order==5) cntr(MRs_order==5);...
             cntr(HOCs_order>=6) cntr(Tfactories_order>=6) cntr(gt1genes_order>=6) cntr(gt2genes_order>=6) cntr(commonTFs_order>=6) cntr(MRs_order>=6)];
         
HOC_table_real_expand = [cntr(HOCs_order==3) cntr(Tfactories_order==3) cntr(gt1genes_order==3) cntr(gt2genes_order==3) cntr(commonTFs_order==3) cntr(MRs_order==3);...
             cntr(HOCs_order==4) cntr(Tfactories_order==4) cntr(gt1genes_order==4) cntr(gt2genes_order==4) cntr(commonTFs_order==4) cntr(MRs_order==4);...
             cntr(HOCs_order==5) cntr(Tfactories_order==5) cntr(gt1genes_order==5) cntr(gt2genes_order==5) cntr(commonTFs_order==5) cntr(MRs_order==5);...
             cntr(HOCs_order==6) cntr(Tfactories_order==6) cntr(gt1genes_order==6) cntr(gt2genes_order==6) cntr(commonTFs_order==6) cntr(MRs_order==6);...
             cntr(HOCs_order==7) cntr(Tfactories_order==7) cntr(gt1genes_order==7) cntr(gt2genes_order==7) cntr(commonTFs_order==7) cntr(MRs_order==7);...
             cntr(HOCs_order==8) cntr(Tfactories_order==8) cntr(gt1genes_order==8) cntr(gt2genes_order==8) cntr(commonTFs_order==8) cntr(MRs_order==8);...
             cntr(HOCs_order==9) cntr(Tfactories_order==9) cntr(gt1genes_order==9) cntr(gt2genes_order==9) cntr(commonTFs_order==9) cntr(MRs_order==9);...
             cntr(HOCs_order==10) cntr(Tfactories_order==10) cntr(gt1genes_order==10) cntr(gt2genes_order==10) cntr(commonTFs_order==10) cntr(MRs_order==10);...
             cntr(HOCs_order==11) cntr(Tfactories_order==11) cntr(gt1genes_order==11) cntr(gt2genes_order==11) cntr(commonTFs_order==11) cntr(MRs_order==11)];

[val, idx] = max(cellfun("length", edgeTable_w_genes_TF.chrNum))

for i = 1:val
    if type == 1
        samp_num_fib(i,1) = cntr(Tfactories_order==i);
    elseif type == 2
        samp_num_gm(i,1) = cntr(Tfactories_order==i);
    end
end

         