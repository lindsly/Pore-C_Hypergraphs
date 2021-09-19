%% Load dataset and convert to matrix
load('NlaIII_digest_IRFibroblast_006_reads.mat') % Fibroblast experiment V2 pairs table
binSize = 1e6;
NlaIII_digest_IRFibrobalst_006_reads.posA = ceil(NlaIII_digest_IRFibrobalst_006_reads.posA/binSize);
NlaIII_digest_IRFibrobalst_006_reads.posB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posB/binSize);
NlaIII_digest_IRFibrobalst_006_reads = sortrows(NlaIII_digest_IRFibrobalst_006_reads, 'read_id', 'ascend');
porecMatrix = table2array(NlaIII_digest_IRFibrobalst_006_reads);
% 
% load('public_gm12878_reads.mat')
% binSize = 1e6;
% public_gm12878_reads.posA = ceil(public_gm12878_reads.posA/binSize);
% public_gm12878_reads.posB = ceil(public_gm12878_reads.posB/binSize);
% public_gm12878_strands = sortrows(public_gm12878_reads, 'read_id', 'ascend');
% porecMatrix = table2array(public_gm12878_reads);

%% Add chromosome offsets
chrLength = [249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
    134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58];
for chr = 2:24
    binAdd = sum(chrLength(1:chr-1));
    binIdxL = porecMatrix(:, 2)==chr;
    porecMatrix(binIdxL, 3) = porecMatrix(binIdxL, 3)+binAdd; 
    binIdxR = porecMatrix(:, 4)==chr;
    porecMatrix(binIdxR, 5) = porecMatrix(binIdxR, 5)+binAdd;
end

%% Threshold pairwise contacts removing noise
[~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 3 5]), 'rows');
porecMatrix = porecMatrix(uniqReadContacts, :);
pairedContacts = porecMatrix(:, [3 5]);
[~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
contactHist = accumarray(contactIdx, 1);
contactCount = contactHist(contactIdx);
eps = prctile(contactCount, 85); % 85th percentile 
porecMatrix = porecMatrix(contactCount>=eps, :); % SAVE THIS for future runs

i_chr = 23;
%% Extract a chromosome of interest
chr = i_chr; % Choose chromosome 22
chrIdx = logical((porecMatrix(:, 2)==chr).*(porecMatrix(:, 4)==chr));
porecMatrixChr = porecMatrix(chrIdx, :);
porecMatrixChr(:, [3 5]) = porecMatrixChr(:, [3 5])-sum(chrLength(1:chr-1));

%% Create incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrixChr(:, 1)); % Reindex the ReadID
porecMatrixChr(:, 1) = readIDReIdx;
incidenceMatrixA = sparse(porecMatrixChr(:, 3), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrixB = sparse(porecMatrixChr(:, 5), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = unique(incidenceMatrix', 'rows')';

%% Compute order frequency
edgeOrders = sum(incidenceMatrix, 1);
maxOrder = max(edgeOrders);
edgeOrderFrequency = zeros(maxOrder, 1);
for i = 1:maxOrder
    edgeOrderFrequency(i) = nnz(edgeOrders==i);
end

%% Extract all higher-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3

%% Create CSV file for PAOHvis
CSVIncidenceMatrix = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
    end 
end
label = i_chr*ones(size(CSVIncidenceMatrix, 1), 1);
        
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%02d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label);
writetable(CSVIncidenceMatrixTable, ['Intra_Chromosome_',num2str(i_chr),'_GM.csv'])
