%% Fibroblast
load('edgeTables_for_summaryTable.mat')

% Total HOC
orderCount1_Fib = orderCount(4, edgeTable_Fibroblast.chrNum);

% Transcption factories candidates
orderCount2_Fib = orderCount(4, edgeTable_Fibroblast_w_genes_TF.chrNum);

% gene>=1 
num1GenesGT1 = cellfun(@(x) sum(x>=1)>=1, edgeTable_Fibroblast_w_genes_TF.geneExpress);
edgeTable_Fibroblast_w_genes_TF_num1 = edgeTable_Fibroblast_w_genes_TF(num1GenesGT1, :);
orderCount3_Fib = orderCount(4, edgeTable_Fibroblast_w_genes_TF_num1.chrNum);

% gene>=2
num2GenesGT1 = cellfun(@(x) sum(x>=1)>=2, edgeTable_Fibroblast_w_genes_TF.geneExpress);
edgeTable_Fibroblast_w_genes_TF_num2 = edgeTable_Fibroblast_w_genes_TF(num2GenesGT1, :);
orderCount4_Fib = orderCount(4, edgeTable_Fibroblast_w_genes_TF_num2.chrNum);

% common TF
numGenesComTF = cellfun(@(x) size(x,1)>0, edgeTable_Fibroblast_w_genes_TF_num2.comTF);
edgeTable_Fibroblast_w_genes_TF_num2_comTF = edgeTable_Fibroblast_w_genes_TF_num2(numGenesComTF, :);
orderCount5_Fib = orderCount(4, edgeTable_Fibroblast_w_genes_TF_num2_comTF.chrNum);

%% GM12878
% Total HOC
orderCount1_GM = orderCount(4, edgeTable_GM12878.chrNum);

% Transcption factories candidates
orderCount2_GM = orderCount(4, edgeTable_GM12878_w_genes_TF.chrNum);

% gene>=1 
num1GenesGT1 = cellfun(@(x) sum(x>=1)>=1, edgeTable_GM12878_w_genes_TF.geneExpress);
edgeTable_GM12878_w_genes_TF_num1 = edgeTable_GM12878_w_genes_TF(num1GenesGT1, :);
orderCount3_GM = orderCount(4, edgeTable_GM12878_w_genes_TF_num1.chrNum);

% gene>=2
num2GenesGT1 = cellfun(@(x) sum(x>=1)>=2, edgeTable_GM12878_w_genes_TF.geneExpress);
edgeTable_GM12878_w_genes_TF_num2 = edgeTable_GM12878_w_genes_TF(num2GenesGT1, :);
orderCount4_GM = orderCount(4, edgeTable_GM12878_w_genes_TF_num2.chrNum);

% common TF
numGenesComTF = cellfun(@(x) size(x,1)>0, edgeTable_GM12878_w_genes_TF_num2.comTF);
edgeTable_GM12878_w_genes_TF_num2_comTF = edgeTable_GM12878_w_genes_TF_num2(numGenesComTF, :);
orderCount5_GM = orderCount(4, edgeTable_GM12878_w_genes_TF_num2_comTF.chrNum);

%% Functions
function count = orderCount(size, tableColumn)
    count = zeros(size, 1);
    for i = 1:size-1
        count(i) = nnz(cellfun('length', tableColumn)==i+2);
    end
    count(end) = nnz(cellfun('length', tableColumn)>=size+2);
end


