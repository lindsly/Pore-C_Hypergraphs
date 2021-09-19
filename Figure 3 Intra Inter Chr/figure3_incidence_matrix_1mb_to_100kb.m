%% Load dataset and convert to matrix
load('NlaIII_digest_IRFibroblast_006_reads.mat')  % Fibroblast experiment V2 pairs table
binSize = 1e5;
binSize1MB = 1e6;
NlaIII_digest_IRFibrobalst_006_reads.posA1MB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posA/binSize1MB);
NlaIII_digest_IRFibrobalst_006_reads.posB1MB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posB/binSize1MB);
NlaIII_digest_IRFibrobalst_006_reads.posA = ceil(NlaIII_digest_IRFibrobalst_006_reads.posA/binSize);
NlaIII_digest_IRFibrobalst_006_reads.posB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posB/binSize);
NlaIII_digest_IRFibrobalst_006_reads = sortrows(NlaIII_digest_IRFibrobalst_006_reads, 'read_id', 'ascend');
porecMatrix = table2array(NlaIII_digest_IRFibrobalst_006_reads);

%% Add chromosome offsets
chrLength = [2490 2422 1983 1903 1816 1709 1594 1452 1384 1338 1351 1333 ...
    1144 1071 1020 904 833 804 587 645 468 509 1561 573];
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
porecMatrix = porecMatrix(contactCount>=eps, :);

%% Extract a chromosome of interest
chr = 22; % Choose chromosome 22
chrIdx = logical((porecMatrix(:, 2)==chr).*(porecMatrix(:, 4)==chr));
porecMatrixChr = porecMatrix(chrIdx, :);
porecMatrixChr(:, [3 5]) = porecMatrixChr(:, [3 5])-sum(chrLength(1:chr-1));

%% Extract 1mb higher-order contacts
highContact = [19 21 22];
idx1 = ismember(porecMatrixChr(:, 6), highContact);
idx2 = ismember(porecMatrixChr(:, 7), highContact);
idx = logical(idx1.*idx2);
porecMatrixChr = porecMatrixChr(idx, :);

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrixChr(:, 1)); % Reindex the ReadID
porecMatrixChr(:, 1) = readIDReIdx;
incidenceMatrixA = sparse(porecMatrixChr(:, 3), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrixB = sparse(porecMatrixChr(:, 5), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = unique(incidenceMatrix', 'rows')';

%% Extract all higer-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3

% Selected for having loci in all three 1 Mb bins
idx_ex = [1 5 10 11 12 14];
highOrderContacts = highOrderContacts(idx_ex);

%% Create CSV file for PAOHvis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        if contact(j) >= 181 && contact(j) <= 190
            label{end+1} = 'Locus L19'; %#ok<*SAGROW>
        elseif contact(j) >= 201 && contact(j) <= 210
            label{end+1} = 'Locus L21';
        else 
            label{end+1} = 'Locus L22';
        end
    end 
end 
        
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%03d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label');
writetable(CSVIncidenceMatrixTable, 'L19_L21_L22_100kb_zoom.csv')

