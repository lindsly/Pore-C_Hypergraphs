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
