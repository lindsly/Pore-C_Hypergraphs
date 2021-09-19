% Load aligned pore-c data

% type 1 is IMR-90 and IR fibroblast, type 2 is GM12878
type = 1;

if type == 1 && exist('aligns_1234','var') == 0
    load('aligns_1234.mat')
    full_align_table = aligns_1234;
    
    full_align_table = sortrows(full_align_table, 'read_idx', 'ascend');
    % Remove the experiment number (already factored into the read_idx column)
    alignPorecMatrix = table2array(full_align_table(:,2:end)); 
elseif type == 2 && exist('public_gm12878_full_alignment_table','var') == 0
    load('public_gm12878_full_alignment_table.mat')
    full_align_table = public_gm12878_full_alignment_table;
    full_align_table.strand(ismember(full_align_table.strand,'+'))='1';
    full_align_table.strand(ismember(full_align_table.strand,'-'))='0';
    full_align_table.strand = double(full_align_table.strand);
    full_align_table.strand(ismember(full_align_table.strand,3))=1;
    full_align_table.strand(ismember(full_align_table.strand,4))=0;
    
    full_align_table = sortrows(full_align_table, 'read_idx', 'ascend');
    alignPorecMatrix = table2array(full_align_table(:, [3 6 7 8 9 14 20 21]));
end
% full_align_table = sortrows(full_align_table, 'read_idx', 'ascend');
% alignPorecMatrix = table2array(full_align_table(:, [3 6 7 8 9 14 20 21]));
alignPorecMatrix = alignPorecMatrix(alignPorecMatrix(:, 6)>=40, :);
alignPorecMatrix = [alignPorecMatrix, alignPorecMatrix(:, 3), alignPorecMatrix(:, 4)];

% Load gene table data
geneLoc = readtable('biomart_hg38_protein_coding.txt', ...
    'ReadVariableNames', true, 'PreserveVariableNames', true);
geneMatrix = table2array(geneLoc(:, [2 3 4]));
allGenes = table2cell(geneLoc(:, 1));

% Load atac-seq data
if type == 1 && exist('imr90_atac','var') == 0
    load('GM12878_IMR90_ATAC.mat')
    DNAse = imr90_atac;
    clear imr90_atac gm12878_atac
elseif type == 2 && exist('gm12878_atac','var') == 0
    load('GM12878_IMR90_ATAC.mat')
    DNAse = gm12878_atac;
    clear imr90_atac gm12878_atac
end
DNAMatrix = table2array(DNAse(:, [1 2 3 7]));
signal = DNAse.signalValue;

% Load chip-seq data
if type == 1 && exist('imr90_atac','var') == 0
    load('GM12878_IMR90_ChIP_Pol2.mat')
    chipSeq = imr90_pol2;
    clear imr90_pol2 gm12878_pol2
elseif type == 2 && exist('gm12878_atac','var') == 0
    load('GM12878_IMR90_ChIP_Pol2.mat')
    chipSeq = gm12878_pol2;
    clear imr90_pol2 gm12878_pol2
end
chipSeqMatrix = table2array(chipSeq(:, [1 2 3 7]));
signal2 = chipSeq.signalValue;

% Load enhancer data
if type == 1 && exist('imr90_atac','var') == 0
    load('GM12878_IMR90_enhancers.mat')
    enhancer = imr90_enhancer;
    clear imr90_enhancer gm12878_enhancer
elseif type == 2 && exist('gm12878_atac','var') == 0
    load('GM12878_IMR90_enhancers.mat')
    enhancer = gm12878_enhancer;
    clear imr90_enhancer gm12878_enhancer
end
enhancerMatrix = table2array(enhancer);
signal3 = enhancer.signalValue;

% Add chr offset
chrLength = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, ...
    159345973, 145138636, 138394717, 133797422, 135086622, ...
    133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, ...
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415];

for chr = 2:24
    binAdd = sum(chrLength(1:chr-1));
    alignBinIdx = alignPorecMatrix(:, 2)==chr;
    alignPorecMatrix(alignBinIdx, 3) = alignPorecMatrix(alignBinIdx, 3)+binAdd;
    alignPorecMatrix(alignBinIdx, 4) = alignPorecMatrix(alignBinIdx, 4)+binAdd;
    
    geneBinIdx = geneMatrix(:, 1)==chr;
    geneMatrix(geneBinIdx, [2 3]) = geneMatrix(geneBinIdx, [2 3])+binAdd;
    
    DNABinIdx = DNAMatrix(:, 1)==chr;
    DNAMatrix(DNABinIdx, [2 3]) = DNAMatrix(DNABinIdx, [2 3])+binAdd;
    
    chipSeqBinIdx = chipSeqMatrix(:, 1)==chr;
    chipSeqMatrix(chipSeqBinIdx, [2 3]) = chipSeqMatrix(chipSeqBinIdx, [2 3])+binAdd;
    
    enhancerBinIdx = enhancerMatrix(:, 1)==chr;
    enhancerMatrix(enhancerBinIdx, [2 3]) = enhancerMatrix(enhancerBinIdx, [2 3])+binAdd;
end

% extract high-order contacts
[~, ~, readIDReIdx] = unique(alignPorecMatrix(:, 1)); % Reindex the ReadID
alignPorecMatrix = [alignPorecMatrix, readIDReIdx]; 
alignPorecMatrix = alignPorecMatrix(:, [11 2 3 4 5 6 7 8 9 10 1]);

%% Random Sampling
% 3rd order
num_3rd = 379165*3;
idx_sample = datasample(1:size(alignPorecMatrix,1),num_3rd,'Replace',false);
alignPorecMatrix = alignPorecMatrix(idx_sample,:);

synth_idx = (1:379165)'; %
synth_idx = repmat(synth_idx,1,3)';
synth_idx = synth_idx(:);
alignPorecMatrix(:,1) = synth_idx;
alignPorecMatrix(:,11) = synth_idx;


%% Continue with normal code
contactHist = accumarray(alignPorecMatrix(:, 1), 1);
contactCount = contactHist(alignPorecMatrix(:, 1));
alignPorecMatrix = alignPorecMatrix(contactCount>=3, :);
[~, ~, readIDReIdx] = unique(alignPorecMatrix(:, 1)); % Reindex the ReadID
alignPorecMatrix(:, 1) = readIDReIdx;
alignPorecMatrix(:, 3) = alignPorecMatrix(:, 3)-5000;
alignPorecMatrix(:, 4) = alignPorecMatrix(:, 4)+5000;

% atac-seq, chip-seq, gene
n = readIDReIdx(end);
readID = zeros(n, 1);
origStartEnd = cell(n, 1);
chrNum = cell(n, 1);
startEndTable = cell(n, 1);
strandTable = cell(n, 1);
fragmentStartEnd = cell(n,1);
mappingQualityTable = cell(n, 1);

geneEdges = cell(n, 1);
DNASignal = cell(n, 1);
chipSeqSignal = cell(n, 1);
enhancerSignal = cell(n, 1);


tic
parfor i = 1:n
    readIDIdx = alignPorecMatrix(:, 1)==i;
    startEnd = alignPorecMatrix(readIDIdx, [3 4]);
    strand = alignPorecMatrix(readIDIdx, 5);
    mappingQuality = alignPorecMatrix(readIDIdx, 6);
    
    genes = cell(0, size(startEnd, 1));
    dnaSignal = zeros(1, size(startEnd, 1));
    chipSignal = zeros(1, size(startEnd, 1));
    enSignal = zeros(1, size(startEnd, 1));
    
    for j = 1:size(startEnd, 1)
        geneIdx = ~(geneMatrix(:, 2)>startEnd(j, 2)|geneMatrix(:, 3)<startEnd(j, 1))...
            & ((startEnd(j, 1)<=geneMatrix(:, 2)&startEnd(j, 2)<=geneMatrix(:, 3)|...
            (startEnd(j, 1)>=geneMatrix(:, 2)&startEnd(j, 2)>=geneMatrix(:, 3)))|...
            (startEnd(j, 1)<=geneMatrix(:, 2)&startEnd(j, 2)>=geneMatrix(:, 3))|...
            (startEnd(j, 1)>=geneMatrix(:, 2)&startEnd(j, 2)<=geneMatrix(:, 3)));
        DNAIdx = ~(DNAMatrix(:, 2)>startEnd(j, 2)|DNAMatrix(:, 3)<startEnd(j, 1))...
            & ((startEnd(j, 1)<=DNAMatrix(:, 2)&startEnd(j, 2)<=DNAMatrix(:, 3)|...
            (startEnd(j, 1)>=DNAMatrix(:, 2)&startEnd(j, 2)>=DNAMatrix(:, 3)))|...
            (startEnd(j, 1)<=DNAMatrix(:, 2)&startEnd(j, 2)>=DNAMatrix(:, 3))|...
            (startEnd(j, 1)>=DNAMatrix(:, 2)&startEnd(j, 2)<=DNAMatrix(:, 3)));
        chiSeqIdx = ~(chipSeqMatrix(:, 2)>startEnd(j, 2)|chipSeqMatrix(:, 3)<startEnd(j, 1))...
            & ((startEnd(j, 1)<=chipSeqMatrix(:, 2)&startEnd(j, 2)<=chipSeqMatrix(:, 3)|...
            (startEnd(j, 1)>=chipSeqMatrix(:, 2)&startEnd(j, 2)>=chipSeqMatrix(:, 3)))|...
            (startEnd(j, 1)<=chipSeqMatrix(:, 2)&startEnd(j, 2)>=chipSeqMatrix(:, 3))|...
            (startEnd(j, 1)>=chipSeqMatrix(:, 2)&startEnd(j, 2)<=chipSeqMatrix(:, 3)));
        enhancerSeqIdx = ~(enhancerMatrix(:, 2)>startEnd(j, 2)|enhancerMatrix(:, 3)<startEnd(j, 1))...
            & ((startEnd(j, 1)<=enhancerMatrix(:, 2)&startEnd(j, 2)<=enhancerMatrix(:, 3)|...
            (startEnd(j, 1)>=enhancerMatrix(:, 2)&startEnd(j, 2)>=enhancerMatrix(:, 3)))|...
            (startEnd(j, 1)<=enhancerMatrix(:, 2)&startEnd(j, 2)>=enhancerMatrix(:, 3))|...
            (startEnd(j, 1)>=enhancerMatrix(:, 2)&startEnd(j, 2)<=enhancerMatrix(:, 3)));
        if nnz(geneIdx) > 0 
            loc = find(geneIdx>0);
            for jj = 1:nnz(geneIdx)
                genes{end+1, j} = allGenes{loc(jj)}; %#ok<SAGROW>
            end
        end
        if nnz(DNAIdx) > 0 
            dnaSignal(j) = sum(signal(DNAIdx));
        end
        if nnz(chiSeqIdx) > 0
            chipSignal(j) = sum(signal2(chiSeqIdx));
        end
        if nnz(enhancerSeqIdx) > 0
            enSignal(j) = sum(signal3(enhancerSeqIdx));
        end
    end 
    readID(i) = unique(alignPorecMatrix(readIDIdx, 11));
    fragmentStartEnd{i} = alignPorecMatrix(readIDIdx,[9 10]);
    origStartEnd{i} = alignPorecMatrix(readIDIdx, [7 8]);
    chrNum{i} = alignPorecMatrix(readIDIdx, 2);
    startEndTable{i} = startEnd;
    strandTable{i} = strand;
    mappingQualityTable{i} = mappingQuality;
    geneEdges{i} = genes;
    DNASignal{i} = dnaSignal;
    chipSeqSignal{i} = chipSignal;
    enhancerSignal{i} = enSignal;
end
toc

edgeTable = table(readID, chrNum, startEndTable, origStartEnd, fragmentStartEnd, strandTable, mappingQualityTable, geneEdges, DNASignal, chipSeqSignal,enhancerSignal);
edgeTable.Properties.VariableNames = {'readID', 'chrNum', 'startEnd', 'origStartEnd', 'fragmentStartEnd', 'strand', 'mappingQuality', 'genes', 'atac-seq', 'chip-seq','enhancer'};

%% %%%%%%%%%% porecAnalysisTransFac3 %%%%%%%%%%
% Load datasets
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
    load('edgeTable_Fibroblast_Combined_enhancer.mat')
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
enhancerIdx = cellfun(@(x) sum(x)>0, edgeTable_w_genes_TF.('enhancer'));
edgeTable_w_genes_TF = edgeTable_w_genes_TF(enhancerIdx, :);
disp(['HOCs with >= 1 enhancer: ', num2str(size(edgeTable_w_genes_TF,1))])

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


HOC_table = [cntr(HOCs_order==3) cntr(Tfactories_order==3) cntr(gt1genes_order==3) cntr(gt2genes_order==3) cntr(commonTFs_order==3) cntr(MRs_order==3);...
             cntr(HOCs_order==4) cntr(Tfactories_order==4) cntr(gt1genes_order==4) cntr(gt2genes_order==4) cntr(commonTFs_order==4) cntr(MRs_order==4);...
             cntr(HOCs_order==5) cntr(Tfactories_order==5) cntr(gt1genes_order==5) cntr(gt2genes_order==5) cntr(commonTFs_order==5) cntr(MRs_order==5);...
             cntr(HOCs_order>=6) cntr(Tfactories_order>=6) cntr(gt1genes_order>=6) cntr(gt2genes_order>=6) cntr(commonTFs_order>=6) cntr(MRs_order>=6)];


